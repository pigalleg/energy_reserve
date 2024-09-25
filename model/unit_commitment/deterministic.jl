using JuMP
using Gurobi
# include("../utils.jl")
# include("../post_processing.jl")
include("./utils.jl")

function DUC(gen_df, loads, gen_variable, mip_gap)
    # model = direct_model(Gurobi.Optimizer(GRB_ENV ))
    
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", mip_gap)
    # set_optimizer_attribute(model, "LogFile", "./output/log_file.txt")
    set_optimizer_attribute(model, "OutputFlag", 0)
    # model = Model(HiGHS.Optimizer)
    # set_optimizer_attribute(model, "mip_rel_gap", mip_gap)
    sets = get_sets(gen_df, loads)
    G = sets.G
    G_thermal = sets.G_thermal
    G_var = sets.G_var
    G_nonvar = sets.G_nonvar
    G_nt_nonvar = sets.G_nt_nonvar
    T = sets.T
    T_red = sets.T_red
    # G, G_thermal, _, G_var, G_nonvar, G_nt_nonvar = create_generators_sets(gen_df)
    # T, T_red = create_time_sets(loads)
    @variables(model, begin
        GEN[G, T]  >= 0     # generation
        COMMIT[G_thermal, T], Bin # commitment status (Bin=binary)
        START[G_thermal, T], Bin  # startup decision
        SHUT[G_thermal, T], Bin   # shutdown decision
    end)
              
  # Objective function
      # Sum of variable costs + start-up costs for all generators and time periods
      # TODO: add delta_T
    @expression(model, StartCost,
        sum(gen_df[gen_df.r_id .== g,:start_cost_per_mw][1]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*START[g,t] for g in G_thermal for t in T)
    )

    @expression(model, OperationalCost,
        sum((gen_df[gen_df.r_id .== g,:heat_rate_mmbtu_per_mwh][1]*gen_df[gen_df.r_id .== g,:fuel_cost][1] + gen_df[gen_df.r_id .== g,:var_om_cost_per_mwh][1])*GEN[g,t] for g in G_nonvar for t in T) +
        sum(gen_df[gen_df.r_id .== g,:var_om_cost_per_mwh][1]*GEN[g,t]  for g in G_var for t in T) + 
        sum(gen_df[gen_df.r_id .== g,:fixed_om_cost_per_mw_per_hour][1]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*COMMIT[g,t] for g in G_thermal for t in T) + 
        sum(gen_df[gen_df.r_id .== g,:fixed_om_cost_per_mw_per_hour][1]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1] for g in G_nt_nonvar for t in T)
    )
    @expression(model, OPEX,
        model[:StartCost] + model[:OperationalCost]
    )

    @objective(model, Min,
        model[:OPEX]
    )
    # Demand balance constraint (supply must = demand in all time periods)
    # Expression is constructed to reuse during ED
    @expression(model, SupplyDemand[t in T],
        sum(GEN[g,t] for g in G)
    )
    @constraint(model, SupplyDemandBalance[t in T], 
        SupplyDemand[t] == loads[loads.hour .== t,:demand][1]
    )

    # Capacity constraints 
    # 1. thermal generators requiring commitment
    @constraint(model, Cap_thermal_min[g in G_thermal, t in T], 
        GEN[g,t] >= COMMIT[g, t]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:min_power][1] 
    ) 
    @constraint(model, Cap_thermal_max[g in G_thermal, t in T], 
        GEN[g,t] <= COMMIT[g, t]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1]
    ) 

    # 2. non-variable generation not requiring commitment
    @constraint(model, Cap_nt_nonvar[g in G_nt_nonvar, t in T], 
        GEN[g,t] <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]
    )

    # 3. variable generation, accounting for hourly capacity factor
    # TODO: The way this constraint is declared does not follow general style
    # Needs to be redefined at each ED
    @constraint(model, Cap_var[g in 1:nrow(gen_variable)], 
            GEN[gen_variable[g,:r_id], gen_variable[g,:hour] ] <= 
                        gen_variable[g,:cf] *
                        gen_variable[g,:existing_cap_mw]
                    )

    # Unit commitment constraints
    # 1. Minimum up time
    @constraint(model, Startup[g in G_thermal, t in T],
        COMMIT[g, t] >= sum(START[g, tt] for tt in intersect(T, (t-gen_df[gen_df.r_id .== g,:up_time][1]):t))
    )

    # 2. Minimum down time
    @constraint(model, Shutdown[g in G_thermal, t in T],
        1-COMMIT[g, t] >= sum(SHUT[g, tt] for tt in intersect(T, (t-gen_df[gen_df.r_id .== g,:down_time][1]):t))
    )

    # 3. Start up/down logic
    @constraint(model, CommitmentStatus[g in G_thermal, t in T_red],
        COMMIT[g,t+1] - COMMIT[g,t] == START[g,t+1] - SHUT[g,t+1]
    )
    return model
end

function add_storage(model, storage, loads, gen_df, sets)
    T = sets.T 
    T_incr = copy(T)
    pushfirst!(T_incr, T_incr[1]-1)
    S = create_storage_sets(storage)
    
    GEN = model[:GEN]
    # START = model[:START]
    @variables(model, begin
        CH[S,T] >= 0
        DIS[S,T] >= 0
        SOE[S,T_incr] >= 0 # T_incr captures SOE at t = T[1]-1
        M[S,T], Bin
    end)

    # Redefinition of objecive function
    @expression(model, StorageOperationalCost,
        sum(storage[storage.r_id .== s,:var_om_cost_per_mwh][1]*(CH[s,t] + DIS[s,t]) for s in S, t in T)
    )

    OPEX = model[:OPEX]
    remove_variable_constraint(model, :OPEX, false)
    @expression(model, OPEX,
        OPEX + model[:StorageOperationalCost]
    )
    @objective(model, Min,
        model[:OPEX]
    )

    # Redefinition of supply-demand balance expression and constraint
    SupplyDemand = model[:SupplyDemand]
    unregister(model, :SupplyDemand)
    @expression(model, SupplyDemand[t in T],
        SupplyDemand[t] - sum(CH[s,t] - DIS[s,t] for s in S)
    )
    SupplyDemandBalance = model[:SupplyDemandBalance]
    delete.(model, SupplyDemandBalance) # Constraints must be deleted also
    unregister(model, :SupplyDemandBalance)
    @constraint(model, SupplyDemandBalance[t in T], 
        SupplyDemand[t] == loads[loads.hour .== t,:demand][1]
    )

    # Charging-discharging logic
    @constraint(model, ChargeLogic[s in S, t in T],
        CH[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1]*M[s,t]
    )
    @constraint(model, DischargeLogic[s in S, t in T],
        DIS[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1]*(1-M[s,t])
    )
    
    # Storage constraints
    @constraint(model, SOEEvol[s in S, t in T], 
        SOE[s,t] == SOE[s,t-1] + CH[s,t]*storage[storage.r_id .== s,:charge_efficiency][1] - DIS[s,t]/storage[storage.r_id .== s,:discharge_efficiency][1]
    ) #TODO: add delta_T

    @constraint(model, SOEMax[s in S, t in T],
        SOE[s,t] <= storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    @constraint(model, SOEMin[s in S, t in T],
        SOE[s,t] >= storage[storage.r_id .== s,:min_energy_mwh][1]
    )
    @constraint(model, CHMin[s in S, t in T],
        CH[s,t] >= storage[storage.r_id .== s,:existing_cap_mw][1]*storage[storage.r_id .== s,:min_power][1] #TODO: calculation should be done at input file
    )
    @constraint(model, DISMin[s in S, t in T],
        DIS[s,t] >= storage[storage.r_id .== s,:existing_cap_mw][1]*storage[storage.r_id .== s,:min_power][1] #TODO: calculation should be done at input file
    )
    # SOE_T_initial = SOE_0
    @constraint(model, SOEO[s in S], #TODO: replace by T
        SOE[s,T_incr[1]] == storage[storage.r_id .== s,:initial_energy_proportion][1]*storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    @constraint(model, SOEFinal[s in S],
        SOE[s,T[end]] == storage[storage.r_id .== s,:initial_energy_proportion][1]*storage[storage.r_id .== s,:max_energy_mwh][1]
    )
end

function add_ramp_constraints(model, gen_df, sets)
    G = sets.G
    G_thermal = sets.G_thermal
    G_nonthermal = sets.G_nonthermal
    T = sets.T
    T_red = sets.T_red

    # G, G_thermal, G_nonthermal, __, ___, ____ = create_generators_sets(gen_df)

    GEN = model[:GEN]
    COMMIT = model[:COMMIT]
    # T = axes(GEN)[2]
    # T_red = T[1:end-1]
    # New auxiliary variable GENAUX for generation above the minimum output level
    @variable(model, GENAUX[G_thermal, T] >= 0)
    
    # for committed thermal units (only created for thermal generators)
    @constraint(model, AuxGen[g in G_thermal, t in T],
        GENAUX[g,t] == GEN[g,t] - COMMIT[g,t]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:min_power][1]
    )
    
    # Ramp equations for thermal generators (constraining GENAUX)
    @constraint(model, RampUp_thermal[g in G_thermal, t in T_red], 
        GENAUX[g,t+1] - GENAUX[g,t] <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_up_percentage][1]
    )

    @constraint(model, RampDn_thermal[g in G_thermal, t in T_red], 
        GENAUX[g,t] - GENAUX[g,t+1] <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_dn_percentage][1]
    )

    # Ramp equations for non-thermal generators (constraining total generation GEN)
    @constraint(model, RampUp_nonthermal[g in G_nonthermal, t in T_red], 
        GEN[g,t+1] - GEN[g,t] <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_up_percentage][1]
    )

    # @constraint(model, RampDn[i in G, t in T_red], 
    #     GEN[i,t] - GEN[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
    #                              gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1])

    @constraint(model, RampDn_nonthermal[g in G_nonthermal, t in T_red], 
        GEN[g,t] - GEN[g,t+1] <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_dn_percentage][1]
    )
end

function add_storage_reserve_power_constraints(model, storage, sets)
    T = sets.T
    S = create_storage_sets(storage)
    CH = model[:CH]
    DIS = model[:DIS]
    @variables(model, begin
        RESUPDIS[S, T] >= 0
        RESUPCH[S, T] >= 0
        RESDNCH[S, T] >= 0
        RESDNDIS[S, T] >= 0
        U[S,T], Bin
        D[S,T], Bin
    end)

    # Reserve up logic - power constraints
    @constraint(model, ResUpStorageDisCapacityMax[s in S, t in T],
        RESUPDIS[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1] - DIS[s,t]
    )
    @constraint(model, ResUpStorageDisLogic[s in S, t in T], #comment this block to obtain Brunninx
        RESUPDIS[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1]*U[s,t]
    )
    @constraint(model, ResUpStorageChCapacityMax[s in S, t in T],
        RESUPCH[s,t] <= + CH[s,t]
    )
    @constraint(model, ResUpStorageChLogic[s in S, t in T], #comment this block to obtain Brunninx
        CH[s,t] - RESUPCH[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1]*(1-U[s,t])
    )
    # Reserve down logic - power constraints
    @constraint(model, ResDownStorageChCapacityMax[s in S, t in T],
        RESDNCH[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1] - CH[s,t]
    )
    @constraint(model, ResDownStorageChLogic[s in S, t in T], #comment this block to obtain Brunninx
        RESDNCH[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1]*D[s,t]
    )
    @constraint(model, ResDownStorageDisCapacityMax[s in S, t in T],
        RESDNDIS[s,t] <= + DIS[s,t]
    )
    @constraint(model, ResDownStorageDisLogic[s in S, t in T], #comment this block to obtain Brunninx
        DIS[s,t] - RESDNDIS[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1]*(1-D[s,t])
    )
end

function add_reserve_constraints(model, reserve, loads, gen_df, storage::Union{DataFrame, Nothing}, bidirectional_storage_reserve::Bool, storage_envelopes::Bool, naive_envelopes::Bool, thermal_reserve::Bool, μ_up::Union{Int64,Float64}, μ_dn::Union{Int64,Float64}, VRESERVE::Union{Int64,Float64}, sets::NamedTuple)
    G_thermal = sets.G_thermal
    T = sets.T
    T_red = sets.T_red
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]
    G_reserve = G_thermal
    if !isnothing(storage)
        S = create_storage_sets(storage)
        G_reserve = union(G_thermal, S)
    end
    @variables(model, begin
        RESUP[G_reserve, T] >= 0
        RESDN[G_reserve, T] >= 0
    end)
    @objective(model, Min, 
        objective_function(model) + VRESERVE*sum(RESUP[g,t] for g in G_reserve, t in T) + VRESERVE*sum(RESDN[g,t] for g in G_reserve, t in T)
    )

    # (1) Reserves limited by committed capacity of generator
    @constraint(model, ResUpThermal[g in G_thermal, t in T],
        RESUP[g,t] <= COMMIT[g,t]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1] - GEN[g,t]
    )
    @constraint(model, ResDnThermal[g in G_thermal, t in T],
        RESDN[g,t] <= GEN[g,t] - COMMIT[g,t]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:min_power][1]
    ) #TODO: calculation should be done at input file
    
    # (2) Reserves limited by ramp rates #TODO: check if this restrictions make sense
    @constraint(model, ResUpRamp[g in G_thermal, t in T],
        RESUP[g,t] <=  gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_up_percentage][1]
    )
    @constraint(model, ResDnRamp[g in G_thermal, t in T],
        RESDN[g,t] <=  gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_dn_percentage][1]
    )
    # (3) Robust ramp constraints
    @constraint(model, ResUpRampRobust[g in G_thermal, t in T_red],
        GEN[g,t+1] + RESUP[g,t+1] - (GEN[g,t] - RESDN[g,t]) <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_up_percentage][1]
    )
    @constraint(model, ResDnRampRobust[g in G_thermal, t in T_red],
        GEN[g,t] + RESUP[g,t] - (GEN[g,t+1] - RESDN[g,t+1]) <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_dn_percentage][1]
    )
    if !thermal_reserve
        for g in G_thermal, t in T
            fix(RESUP[g,t], 0.0, force = true)
            fix(RESDN[g,t], 0.0, force = true)
        end
    end

    # (3) Storage reserve
    if !isnothing(storage)
        S = create_storage_sets(storage)
        SOE = model[:SOE]
        add_storage_reserve_power_constraints(model, storage, sets)
        RESUPDIS = model[:RESUPDIS]
        RESUPCH = model[:RESUPCH]
        RESDNCH = model[:RESDNCH]
        RESDNDIS = model[:RESDNDIS]

        # Energy constraints # important for μ<1
        @constraint(model, ResUpStorageDisMax[s in S, t in T],
            RESUPDIS[s,t] <= (SOE[s,t]- storage[storage.r_id .== s,:min_energy_mwh][1])*storage[storage.r_id .== s,:discharge_efficiency][1] #TODO: include delta_T
        )
        @constraint(model, ResDownStorageChMax[s in S, t in T],
            RESDNCH[s,t] <= (storage[storage.r_id .== s,:max_energy_mwh][1] - SOE[s,t])/storage[storage.r_id .== s,:charge_efficiency][1] #TODO: include delta_T
        )
        
        @constraint(model, ResUpStorageCapacityMax[s in S, t in T],
            RESUP[s,t] == RESUPDIS[s,t] + RESUPCH[s,t]
        )
        @constraint(model, ResDownStorageCapacityMax[s in S, t in T],
            RESDN[s,t] == RESDNCH[s,t] + RESDNDIS[s,t]
        )

        if !bidirectional_storage_reserve
            for s in S, t in T
                fix(RESUPCH[s,t], 0, force = true)
                fix(RESDNDIS[s,t], 0, force = true)
                fix(U[s,t], 1, force = true)
                fix(D[s,t], 1, force = true)
            end
            remove_variable_constraint(model, :ResUpStorageDisLogic)
            remove_variable_constraint(model, :ResUpStorageChLogic)
            remove_variable_constraint(model, :ResDownStorageChLogic)
            remove_variable_constraint(model, :ResDownStorageDisLogic)
            
        end

        if storage_envelopes # TODO: change this to default when reserves are active
            println("Adding storage envelopes...")
            add_envelope_constraints(model, loads, storage, μ_up, μ_dn, naive_envelopes)
        end
    end
    # (4) Overall reserve requirements
    @constraint(model, ResUpRequirement[t in T],
        sum(RESUP[g,t] for g in G_reserve) >= reserve[reserve.hour .== t,:reserve_up_MW][1]
    )
    @constraint(model, ResDnRequirement[t in T],
        sum(RESDN[g,t] for g in G_reserve) >= reserve[reserve.hour .== t,:reserve_down_MW][1]
    )

end

function add_envelope_constraints(model, loads, storage, μ_up, μ_dn, naive_envelopes = false)
    S = create_storage_sets(storage)
    RESUPCH = model[:RESUPCH]
    RESDNCH = model[:RESDNCH]
    RESUPDIS = model[:RESUPDIS]
    RESDNDIS = model[:RESDNDIS]

    SOE = model[:SOE]
    CH = model[:CH]
    DIS = model[:DIS]
    T, _ =  create_time_sets(loads)
    T_incr = copy(T)
    pushfirst!(T_incr, T_incr[1]-1)
    @variables(model, begin
        SOEUP[S, T_incr] >= 0
        SOEDN[S, T_incr] >= 0
    end)
    if !naive_envelopes
        @constraint(model, SOEUpEvol[s in S, t in T],
            SOEUP[s,t]  == SOEUP[s,t-1] + (CH[s,t] + μ_dn*RESDNCH[s,t])*storage[storage.r_id .== s,:charge_efficiency][1] - (DIS[s,t] - μ_dn*RESDNDIS[s,t])/storage[storage.r_id .== s,:discharge_efficiency][1]
        ) #TODO: add delta_T
        @constraint(model, SOEDnEvol[s in S, t in T], 
            SOEDN[s,t]  == SOEDN[s,t-1] + (CH[s,t] - μ_up*RESUPCH[s,t])*storage[storage.r_id .== s,:charge_efficiency][1] - (DIS[s,t] + μ_up*RESUPDIS[s,t])/storage[storage.r_id .== s,:discharge_efficiency][1]
        ) #TODO: add delta_T
    else
        @constraint(model, SOEUpEvol[s in S, t in T],
            SOEUP[s,t]  == SOE[s,t] + μ_dn*RESDNCH[s,t]*storage[storage.r_id .== s,:charge_efficiency][1] + μ_dn*RESDNDIS[s,t]/storage[storage.r_id .== s,:discharge_efficiency][1]
        ) 
        @constraint(model, SOEDnEvol[s in S, t in T], 
            SOEDN[s,t]  == SOE[s,t] - μ_up*RESUPCH[s,t]*storage[storage.r_id .== s,:charge_efficiency][1] - μ_up*RESUPDIS[s,t]/storage[storage.r_id .== s,:discharge_efficiency][1]
        )
    end
    # SOEUP_T_initial = SOE_T_initial
    @constraint(model, SOEUP_0[s in S],
        SOEUP[s,T_incr[1]] == SOE[s,T_incr[1]]
    )
    # SOEDN_T_initial = SOE_T_initial
    @constraint(model, SOEDN_0[s in S],
        SOEDN[s,T_incr[1]] == SOE[s,T_incr[1]]
    )
    
    # SOEUP, SOEDN <=SOE_max
    @constraint(model, SOEUPMax[s in S, t in T],
        SOEUP[s,t] <= storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    @constraint(model, SOEDNMax[s in S, t in T],
        SOEDN[s,t] <= storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    
    # SOEUP, SOEDN>=SOE_min
    @constraint(model, SOEUPMin[s in S, t in T],
        SOEUP[s,t] >= storage[storage.r_id .== s,:min_energy_mwh][1]
    )
    @constraint(model, SOEDNMin[s in S, t in T],
        SOEDN[s,t] >= storage[storage.r_id .== s,:min_energy_mwh][1]
    )
end

function add_energy_reserve_constraints(model, reserve, loads, gen_df, storage::Union{DataFrame, Nothing}, storage_envelopes::Bool, storage_link_constraint::Bool, thermal_reserve::Bool, μ_up::Union{Int64,Float64}, μ_dn::Union{Int64,Float64}, VRESERVE::Union{Int64,Float64}, sets::NamedTuple)
    #TODO: include diagonal ramp reserves
    G_thermal = sets.G_thermal
    T = sets.T
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]

    G_reserve = G_thermal
    if !isnothing(storage)
        S = create_storage_sets(storage)
        G_reserve = union(G_thermal, S)
        SOE = model[:SOE]
    end

    @variables(model, begin
        ERESUP[G_reserve, j in T, t in T; j <= t] >= 0
        ERESDN[G_reserve, j in T, t in T; j <= t] >= 0
    end)

    @objective(model, Min, 
        objective_function(model) + VRESERVE*sum(ERESUP[g,t,t] for g in G_reserve, t in T) + VRESERVE*sum(ERESDN[g,t,t] for g in G_reserve, t in T)
    )

    # (1) Reserves limited by committed capacity of generator
    @constraint(model, EnergyResUpThermal[g in G_thermal, j in T, t in T; j <= t],
        ERESUP[g, j, t] <= sum(
            COMMIT[g,tt]*gen_df[gen_df.r_id .== g, :existing_cap_mw][1] - GEN[g,tt] for tt in T if (tt >= j)&(tt <= t)
        )
    )
    @constraint(model, EnergyResDownThermal[g in G_thermal, j in T, t in T; j <= t],
        ERESDN[g, j, t] <= sum(
            GEN[g,tt] - COMMIT[g,tt]*gen_df[gen_df.r_id .==g,:existing_cap_mw][1]*gen_df[gen_df.r_id .==g,:min_power][1] for tt in T if (tt >= j)&(tt <= t) #TODO: calculation should be done at input file
        )
    )
    # (2) Reserves limited by ramp rates
    @constraint(model, EnergyResUpRamp[g in G_thermal, j in T, t in T; j <= t],
        ERESUP[g, j, t] <=  sum(
            gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_up_percentage][1] for tt in T if (tt >= j)&(tt <= t) #TODO: calculation should be done at input file
        )
    )
    @constraint(model, EnergyResDnRamp[g in G_thermal, j in T, t in T; j <= t],
        ERESDN[g, j, t] <=  sum(
            gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_dn_percentage][1] for tt in T if (tt >= j)&(tt <= t) #TODO: calculation should be done at input file
        )
    )
    if !thermal_reserve
        for (g,j,t) in [(g,j,t) for g in G_thermal, j in T, t in T if j <= t]
            fix(ERESUP[g,j,t], 0.0, force = true)
            fix(ERESDN[g,j,t], 0.0, force = true)
        end
    end

    # (3) Storage reserve
    if !isnothing(storage)
        # -- begin -- 
        # Following variables are not constrained when energy reserves are used.
        add_storage_reserve_power_constraints(model, storage, sets)
        RESUPDIS = model[:RESUPDIS]
        RESUPCH = model[:RESUPCH]
        RESDNCH = model[:RESDNCH]
        RESDNDIS = model[:RESDNDIS]
        # --end --
        @variables(model, begin
            ERESUPDIS[S, j in T, t in T; j <= t] >= 0
            ERESUPCH[S, j in T, t in T; j <= t] >= 0
            ERESDNCH[S, j in T, t in T; j <= t] >= 0
            ERESDNDIS[S, j in T, t in T; j <= t] >= 0
        end)

        # Power constraints 
        @constraint(model, EnergyResUpStorageDisCapacityMax[s in S, j in T, t in T; j <= t],
            ERESUPDIS[s,j,t] <= sum(
                 RESUPDIS[s,tt] for tt in T if (tt >= j)&(tt <= t)
            )
        )
        @constraint(model, EnergyResUpStorageChCapacityMax[s in S, j in T, t in T; j <= t],
            ERESUPCH[s,j,t] <= sum(
                RESUPCH[s,tt] for tt in T if (tt >= j)&(tt <= t)
            )
        )
        @constraint(model, EnergyResDownStorageChCapacityMax[s in S, j in T, t in T; j <= t],
            ERESDNCH[s,j,t] <= sum(
                RESDNCH[s,tt] for tt in T if (tt >= j)&(tt <= t)
            )
        )
        @constraint(model, EnergyResDownStorageDisCapacityMax[s in S, j in T, t in T; j <= t],
            ERESDNDIS[s,j,t] <= sum(
                RESDNDIS[s,tt] for tt in T if (tt >= j)&(tt <= t)
            )
        )

        # Energy constraints
        @constraint(model, EnergyResUpStorageEnergyMax[s in S, j in T, t in T; j <= t], # important for μ<1
            ERESUPDIS[s,j,t] <= (SOE[s,t]- storage[storage.r_id .== s,:min_energy_mwh][1])*storage[storage.r_id .== s,:discharge_efficiency][1]
        )
        @constraint(model, EnergyResDownStorageEnergyMax[s in S, j in T, t in T; j <= t], # important for μ<1
            ERESDNCH[s,j,t] <= (storage[storage.r_id .== s,:max_energy_mwh][1] - SOE[s,t])/storage[storage.r_id .== s,:charge_efficiency][1]
        )
        # @constraint(model, EnergyResUpStorageDisCapacityMax2[s in S, t in T],
        #     ERESUPDIS[s,t,t] <= 0.5*RESUPDIS[s,t] 
        # )
        # @constraint(model, EnergyResUpStorageChCapacityMax2[s in S, t in T],
        #     ERESUPCH[s,t,t] <= 0.5*RESUPCH[s,t] 
        # )
        # @constraint(model, EnergyResDownStorageChCapacityMax2[s in S, t in T],
        #     ERESDNCH[s,t,t] <= 0.5*RESDNCH[s,t] 
        # )
        # @constraint(model, EnergyResDownStorageDisCapacityMax2[s in S, t in T],
        #     ERESDNDIS[s,t,t] <= 0.5*RESDNDIS[s,t]    
        # )

        # Energy constraints / envelope-like constraints
        @constraint(model, EnergyResUpStorage[s in S, j in T, t in T; j <= t],
            ERESUP[s,j,t] == ERESUPDIS[s,j,t] + ERESUPCH[s,j,t]
        )
        @constraint(model, EnergyResDownStorage[s in S, j in T, t in T; j <= t],
            ERESDN[s,j,t] == ERESDNCH[s,j,t] + ERESDNDIS[s,j,t]
        )

        if storage_envelopes # TODO: change this to default when reserves are active
            println("Adding storage energy envelopes...")
            add_energy_envelope_constraints(model, storage, μ_up, μ_dn, sets)
        end

        if storage_link_constraint
            @constraint(model, EnergyResUpLink[s in S, j in T, t in T; j <= t],
                ERESUPDIS[s,j,t] == sum(ERESUPDIS[s, tt, tt] for tt in T if (tt >= j)&(tt <= t))
            )
            @constraint(model, EnergyResUpLinkBis[s in S, j in T, t in T; j <= t],
                ERESUPCH[s,j,t] == sum(ERESUPCH[s, tt, tt] for tt in T if (tt >= j)&(tt <= t))
            )
            @constraint(model, EnergyResDownLink[s in S, j in T, t in T; j <= t],
                ERESDNCH[s,j,t] == sum(ERESDNCH[s, tt, tt] for tt in T if (tt >= j)&(tt <= t))
            )
            @constraint(model, EnergyResDownLinkBis[s in S, j in T, t in T; j <= t],
                ERESDNDIS[s,j,t] == sum(ERESDNDIS[s, tt, tt] for tt in T if (tt >= j)&(tt <= t))
            )
        end
    end

    # (4) Overall reserve requirements
    @constraint(model, EnergyResUpRequirement[j in T, t in T; j <= t],
        sum(ERESUP[i, j, t] for i in G_reserve) >= reserve[(reserve.i_hour .== j).&(reserve.t_hour .== t),:reserve_up_MW][1]
    )
 
    @constraint(model, EnerResDnRequirement[j in T, t in T; j <= t],
        sum(ERESDN[i, j, t] for i in G_reserve) >= reserve[(reserve.i_hour .== j).&(reserve.t_hour .== t),:reserve_down_MW][1]
    )

end


function add_energy_envelope_constraints(model, storage, μ_up, μ_dn, sets)
    SOE = model[:SOE]
    S = axes(SOE)[1]
    T = sets.T
    ERESUPCH = model[:ERESUPCH]
    ERESDNCH = model[:ERESDNCH]
    ERESUPDIS = model[:ERESUPDIS]
    ERESDNDIS = model[:ERESDNDIS]
    
    @variables(model, begin
        ESOEUP[S, j in T, t in T; j <= t] >= 0
        ESOEDN[S, j in T, t in T; j <= t] >= 0
    end)
    @constraint(model, ESOEUpEvol[s in S, j in T, t in T; j <= t],
        ESOEUP[s,j,t]  == SOE[s,t] + μ_dn*ERESDNCH[s,j,t]*storage[storage.r_id .== s,:charge_efficiency][1] + μ_dn*ERESDNDIS[s,j,t]/storage[storage.r_id .== s,:discharge_efficiency][1]
    )
    @constraint(model, ESOEDnEvol[s in S, j in T, t in T; j <= t], 
        ESOEDN[s,j,t]  == SOE[s,t] - μ_up*ERESUPCH[s,j,t]*storage[storage.r_id .== s,:charge_efficiency][1] - μ_up*ERESUPDIS[s,j,t]/storage[storage.r_id .== s,:discharge_efficiency][1]
    )
    
    # SOEUP, SOEDN <=SOE_max
    @constraint(model, ESOEUPMax[s in S, j in T, t in T; j <= t],
        ESOEUP[s,j,t] <= storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    @constraint(model, ESOEDNMax[s in S, j in T, t in T; j <= t],
        ESOEDN[s,j,t] <= storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    
    # SOEUP, SOEDN>=SOE_min
    @constraint(model, ESOEUPMin[s in S, j in T, t in T; j <= t],
        ESOEUP[s,j,t] >= storage[storage.r_id .== s,:min_energy_mwh][1]
    )
    @constraint(model, ESOEDNMin[s in S, j in T, t in T; j <= t],
        ESOEDN[s,j,t] >= storage[storage.r_id .== s,:min_energy_mwh][1]
    )
end

function add_MGA_constraints(model, reference_solution)
    #deprecated
  
    reference_objective_function = objective_value(model)
    
    mip_gap = get_optimizer_attribute(model,"MIPGap")
    primal = objective_bound(model)
    dual = dual_objective_value(model)

    T = axes(model[:ResUpThermal])[2]
    G_thermal = axes(model[:ResUpThermal])[1]
    S = axes(model[:ResUpStorageCapacityMax])[1]
    RESUP = model[:RESUP]
    RESDN = model[:RESDN]
    COMMIT = model[:COMMIT]
    # LOL  = model[:LOL]
    # remove_variable_constraint(model, :OVMax)
    # remove_variable_constraint(model, :OVMin)
    # remove_variable_constraint(model, :ResUpThermalMin)
    # remove_variable_constraint(model, :ResDownThermalMin)
    # remove_variable_constraint(model, :CommitmentMin)
    # remove_variable_constraint(model, :ResUpStorageMax)
    # remove_variable_constraint(model, :ResUpStorageMin)

    @constraint(model, OVMax,
        objective_function(model) <= (1+mip_gap)*primal
        # objective_function(model) <= primal
    )
    @constraint(model, OVMin,
        objective_function(model) >= (1-mip_gap)*primal
        # objective_function(model) >= dual
    )

    # set_optimizer_attribute(model, "PoolGap", mip_gap)
    
    # @constraint(model, ResUpStorageMin[s in S, t in T],
    #     RESUP[s,t] >=  (filter([:r_id, :hour] => ((x,y) -> (x == s)*(y==t)), reference_solution.reserve).reserve_up_MW)[1]
    # )
    # T_ = 169:176
    # T_ = [172]
    # @constraint(model, ResUpStorageMax[t in T_],
    #     RESUP[101,t] <=  0
    # )
    # @constraint(model, ResUpStorageMax[s in S, t in T],
    #     RESUP[s,t] <=  (filter([:r_id, :hour] => ((x,y) -> (x == s)*(y==t)), reference_solution.reserve).reserve_up_MW)[1]
    # )
    # @constraint(model, ResDownStorageMin[s in S, t in T],,
    #     sum(RESDN[i,t] for i in G_thermal, t in T) >=  sum(filter(:r_id => x -> x in G_thermal, reference_solution.reserve).reserve_down_MW)
    # )

    # @constraint(model, ResUpThermalMin,
    #     sum(RESUP[i,t] for i in G_thermal, t in T) >=  sum(filter(:r_id => x -> x in G_thermal, reference_solution.reserve).reserve_up_MW)
    # )
    # @constraint(model, ResDownThermalMin,
    #     sum(RESDN[i,t] for i in G_thermal, t in T) >=  sum(filter(:r_id => x -> x in G_thermal, reference_solution.reserve).reserve_down_MW)
    # )

    # @constraint(model, ResUpThermalMin[t in T],
    #     sum(RESUP[i,t] for i in G_thermal) >=  sum(filter([:r_id, :hour] => ((x,y) -> (x in G_thermal)*(y==t)), reference_solution.reserve).reserve_up_MW)
    # )
    # @constraint(model, ResDownThermalMin[t in T],
    #     sum(RESDN[i,t] for i in G_thermal) >=  sum(filter([:r_id, :hour] => ((x,y) -> (x in G_thermal)*(y==t)), reference_solution.reserve).reserve_down_MW)
    # )

    
    # @constraint(model, CommitmentMin,
    #     sum(COMMIT[i,t] for i in G_thermal, t in T) >= sum(filter(:r_id => x -> x in G_thermal, reference_solution.generation).commit)
    # )

    # @constraint(model, ResUpThermalMin[i in G_thermal, t in T],
    # RESUP[i,t] >=  sum(filter([:r_id, :hour] => ((x,y) -> (x ==i)*(y==t)), reference_solution.reserve).reserve_up_MW)
    # )
    # @constraint(model, ResDownThermalMin[i in G_thermal, t in T],
    #     RESDN[i,t] >=  sum(filter([:r_id, :hour] => ((x,y) -> (x ==i)*(y==t)), reference_solution.reserve).reserve_down_MW)
    # )
    # @constraint(model, CommitmentMin[i in G_thermal, t in T],
    #     COMMIT[i,t] >= sum(filter([:r_id, :hour] => ((x,y) -> (x ==i)*(y==t)), reference_solution.generation).commit)
    # )
end

function generate_alternative_model(model, reference_solution)
    # Assumes reference_solution is enriched (cf. enrich_dfs())
    println("Adding MGA constraints...")
    # model = JuMP.copy(model)
    # set_optimizer(model, Gurobi.Optimizer) 
    optimize!(model) # Optimizitation is needed since OV will be used as reference.
    # aux = get_solution(uc).scalar.objective_value[1]
    add_MGA_constraints(model, reference_solution)
    # optimize!(model)
    return model
end


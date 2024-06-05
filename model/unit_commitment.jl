using JuMP
using Gurobi
using DataFrames
# using HiGHS
include("./utils.jl")
include("./post_processing.jl")


function remove_variable_constraint(model, key, delete_ = true)
    # Applies for constraints and variables
    println("Removing $key...")
    if !haskey(model, key)
        println("variable not in model")
        return
    end
    if delete_ delete.(model, model[key]) end # Constraints must be deleted also
    unregister(model, key)
end

# TODO change gen_variable => gen_varialbe_df, loads => loads_df
function create_generators_sets(gen_df)
    # Thermal resources for which unit commitment constraints apply
    G_thermal = gen_df[gen_df[!,:up_time] .> 0,:r_id] 
        
    # Non-thermal resources for which unit commitment constraints do NOT apply 
    G_nonthermal = gen_df[gen_df[!,:up_time] .== 0,:r_id]
    
    # Variable renewable resources
    G_var = gen_df[gen_df[!,:is_variable] .== 1,:r_id]
    
    # Non-variable (dispatchable) resources
    G_nonvar = gen_df[gen_df[!,:is_variable] .== 0,:r_id]
    
    # Non-variable and non-thermal resources
    G_nt_nonvar = intersect(G_nonvar, G_nonthermal)
    # Note that G_nt_var = G_var

    # Set of all generators (above are all subsets of this)
    G = gen_df.r_id

    return G, G_thermal, G_nonthermal, G_var, G_nonvar, G_nt_nonvar
end

function create_time_sets(loads)
    return loads.hour, loads.hour[1:end-1]
end

function create_storage_sets(storage)
    return storage.r_id
end
function unit_commitment(gen_df, loads, gen_variable, mip_gap)
    # model = direct_model(Gurobi.Optimizer(GRB_ENV ))
    
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", mip_gap)
    # set_optimizer_attribute(model, "LogFile", "./output/log_file.txt")
    set_optimizer_attribute(model, "OutputFlag", 0)

    # model = Model(HiGHS.Optimizer)
    # set_optimizer_attribute(model, "mip_rel_gap", mip_gap)

    G, G_thermal, _, G_var, G_nonvar, G_nt_nonvar = create_generators_sets(gen_df)
    T, T_red =  create_time_sets(loads)

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
        sum(gen_df[gen_df.r_id .== i,:start_cost_per_mw][1]*gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*START[i,t] for i in G_thermal for t in T)
    )

    @expression(model, OperationalCost,
        sum((gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1]*gen_df[gen_df.r_id .== i,:fuel_cost][1] + gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1])*GEN[i,t] for i in G_nonvar for t in T) +
        sum(gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1]*GEN[i,t]  for i in G_var for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:fixed_om_cost_per_mw_per_hour][1]*gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*COMMIT[i,t] for i in G_thermal for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:fixed_om_cost_per_mw_per_hour][1]*gen_df[gen_df.r_id .== i,:existing_cap_mw][1] for i in G_nt_nonvar for t in T)
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
        sum(GEN[i,t] for i in G)
    )
    @constraint(model, SupplyDemandBalance[t in T], 
        SupplyDemand[t] == loads[loads.hour .== t,:demand][1]
    )

    # Capacity constraints 
    # 1. thermal generators requiring commitment
    @constraint(model, Cap_thermal_min[i in G_thermal, t in T], 
        GEN[i,t] >= COMMIT[i, t]*gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:min_power][1] 
    ) 
    @constraint(model, Cap_thermal_max[i in G_thermal, t in T], 
        GEN[i,t] <= COMMIT[i, t]*gen_df[gen_df.r_id .== i,:existing_cap_mw][1]
    ) 

    # 2. non-variable generation not requiring commitment
    @constraint(model, Cap_nt_nonvar[i in G_nt_nonvar, t in T], 
        GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1]
    )

    # 3. variable generation, accounting for hourly capacity factor
    # TODO: The way this constraint is declared does not follow general style
    # Needs to be redefined at each ED
    @constraint(model, Cap_var[i in 1:nrow(gen_variable)], 
            GEN[gen_variable[i,:r_id], gen_variable[i,:hour] ] <= 
                        gen_variable[i,:cf] *
                        gen_variable[i,:existing_cap_mw]
                    )

    # Unit commitment constraints
    # 1. Minimum up time
    @constraint(model, Startup[i in G_thermal, t in T],
        COMMIT[i, t] >= sum(START[i, tt] for tt in intersect(T, (t-gen_df[gen_df.r_id .== i,:up_time][1]):t))
    )

    # 2. Minimum down time
    @constraint(model, Shutdown[i in G_thermal, t in T],
        1-COMMIT[i, t] >= sum(SHUT[i, tt] for tt in intersect(T, (t-gen_df[gen_df.r_id .== i,:down_time][1]):t))
    )

    # 3. Commitment state
    @constraint(model, CommitmentStatus[i in G_thermal, t in T_red],
        COMMIT[i,t+1] - COMMIT[i,t] == START[i,t+1] - SHUT[i,t+1]
    )
    return model
end

function add_storage(model, storage, loads, gen_df)
    G, G_thermal, __, G_var, G_nonvar, ___ = create_generators_sets(gen_df)
    T, ____ =  create_time_sets(loads)
    T_incr = copy(T)
    pushfirst!(T_incr, T_incr[1]-1)
    S = create_storage_sets(storage)
    big_M = 1000
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
        CH[s,t] <= big_M*M[s,t]
    )
    @constraint(model, DischargeLogic[s in S, t in T],
        DIS[s,t] <= big_M*(1-M[s,t])
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
    @constraint(model, CHMax[s in S, t in T],
        CH[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1]
    )
    @constraint(model, CHMin[s in S, t in T],
        CH[s,t] >= storage[storage.r_id .== s,:existing_cap_mw][1]*storage[storage.r_id .== s,:min_power][1] #TODO: calculation should be done at input file
    )
    @constraint(model, DISMax[s in S, t in T],
        DIS[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1]
    )
    @constraint(model, DISMin[s in S, t in T],
        DIS[s,t] >= storage[storage.r_id .== s,:existing_cap_mw][1]*storage[storage.r_id .== s,:min_power][1] #TODO: calculation should be done at input file
    )
    # SOE_T_initial = SOE_0
    @constraint(model, SOEO[s in S], #TODO: replace by T
        SOE[s,T_incr[1]] == storage[storage.r_id .== s,:initial_energy_proportion][1]*storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    @constraint(model, SOEFinal[s in S],
        SOE[s,T[end]] >= storage[storage.r_id .== s,:initial_energy_proportion][1]*storage[storage.r_id .== s,:max_energy_mwh][1]
    )
end

function add_ramp_constraints(model, loads, gen_df)
    G, G_thermal, G_nonthermal, __, ___, ____ = create_generators_sets(gen_df)
    T, T_red =  create_time_sets(loads)
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]

    # New auxiliary variable GENAUX for generation above the minimum output level
    @variable(model, GENAUX[G_thermal, T] >= 0)
    
    # for committed thermal units (only created for thermal generators)
    @constraint(model, AuxGen[i in G_thermal, t in T],
        GENAUX[i,t] == GEN[i,t] - COMMIT[i, t]*gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:min_power][1]
    )
    
    # Ramp equations for thermal generators (constraining GENAUX)
    @constraint(model, RampUp_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t+1] - GENAUX[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_up_percentage][1]
    )

    @constraint(model, RampDn_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t] - GENAUX[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1]
    )

    # Ramp equations for non-thermal generators (constraining total generation GEN)
    @constraint(model, RampUp_nonthermal[i in G_nonthermal, t in T_red], 
        GEN[i,t+1] - GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_up_percentage][1]
    )

    # @constraint(model, RampDn[i in G, t in T_red], 
    #     GEN[i,t] - GEN[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
    #                              gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1])

    @constraint(model, RampDn_nonthermal[i in G_nonthermal, t in T_red], 
        GEN[i,t] - GEN[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1]
    )
end

function add_reserve_constraints(model, reserve, loads, gen_df, storage::Union{DataFrame, Nothing}, storage_envelopes::Bool, μ_up::Float64, μ_dn::Float64, VRESERVE::Float64)
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]
    _, G_thermal, _, __, ___, ____ = create_generators_sets(gen_df)
    T, T_red =  create_time_sets(loads)
    G_reserve = G_thermal
    if !isnothing(storage)
        S = create_storage_sets(storage)
        G_reserve = union(G_thermal, S)
        SOE = model[:SOE]
        CH = model[:CH]
        DIS = model[:DIS]

    end
    @variables(model, begin
        RESUP[G_reserve, T] >= 0
        RESDN[G_reserve, T] >= 0

    end)
    @objective(model, Min, 
        objective_function(model) + VRESERVE*sum(RESUP[g,t] for g in G_reserve, t in T) + VRESERVE*sum(RESDN[g,t] for g in G_reserve, t in T)
    )

    # (1) Reserves limited by committed capacity of generator
    @constraint(model, ResUpCap[i in G_thermal, t in T],
        RESUP[i,t] <= COMMIT[i,t]*gen_df[gen_df.r_id .==i,:existing_cap_mw][1] - GEN[i,t]
    )
    @constraint(model, ResDnCap[i in G_thermal, t in T],
        RESDN[i,t] <= GEN[i,t] - COMMIT[i,t]*gen_df[gen_df.r_id .==i,:existing_cap_mw][1]*gen_df[gen_df.r_id .==i,:min_power][1]
    ) #TODO: calculation should be done at input file
    
    # (2) Reserves limited by ramp rates
    @constraint(model, ResUpRamp[i in G_thermal, t in T],
        RESUP[i,t] <=  gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_up_percentage][1]
    )
    @constraint(model, ResDnRamp[i in G_thermal, t in T],
        RESDN[i,t] <=  gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1]
    )
    # (3) Robus ramp constraints
    @constraint(model, ResUpRampRobust[i in G_thermal, t in T_red],
        GEN[i,t+1] + RESUP[i,t+1] - (GEN[i,t] - RESDN[i,t]) <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_up_percentage][1]
    )
    @constraint(model, ResDnRampRobust[i in G_thermal, t in T_red],
        GEN[i,t] + RESUP[i,t] - (GEN[i,t+1] - RESDN[i,t+1]) <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1]
    )

    # @constraint(model, ResUpRampZero[i in G_thermal, t in T],
    #     RESUP[i,t] <= 0
    # )
    # @constraint(model, ResDnRampZero[i in G_thermal, t in T],
    #     RESDN[i,t] <= 0
    # )

    # (3) Storage reserve
    if !isnothing(storage)
        @constraint(model, ResUpStorage[s in S, t in T],
            RESUP[s,t] <= (SOE[s,t]- storage[storage.r_id .== s,:min_energy_mwh][1])*storage[storage.r_id .== s,:discharge_efficiency][1] #TODO: include delta_T
        )
        @constraint(model, ResDownStorage[s in S, t in T],
            RESDN[s,t] <= (storage[storage.r_id .== s,:max_energy_mwh][1] - SOE[s,t])/storage[storage.r_id .== s,:charge_efficiency][1] #TODO: include delta_T
        )

        @constraint(model, ResUpStorageCapacityMax[s in S, t in T],
            RESUP[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1] - DIS[s,t] + CH[s,t]
        )
        @constraint(model, ResDownStorageCapacityMax[s in S, t in T],
            RESDN[s,t] <= storage[storage.r_id .== s,:existing_cap_mw][1] - CH[s,t] + DIS[s,t]
        )

        if storage_envelopes
            println("Adding storage envelopes...")
            add_envelope_constraints(model, loads, storage, μ_up, μ_dn)
        end
    end

    # (4) Overall reserve requirements
    @constraint(model, ResUpRequirement[t in T],
        sum(RESUP[i,t] for i in G_reserve) >= reserve[reserve.hour .== t,:reserve_up_MW][1]
    )
    @constraint(model, ResDnRequirement[t in T],
        sum(RESDN[i,t] for i in G_reserve) >= reserve[reserve.hour .== t,:reserve_down_MW][1]
    )
end

function add_envelope_constraints(model, loads, storage, μ_up, μ_dn)
    S = create_storage_sets(storage)
    RESUP = model[:RESUP]
    RESDN = model[:RESDN]
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
    @constraint(model, SOEUpEvol[s in S, t in T], 
        SOEUP[s,t]  == SOEUP[s,t-1] + (CH[s,t] + μ_up*RESDN[s,t])*storage[storage.r_id .== s,:charge_efficiency][1] - DIS[s,t]/storage[storage.r_id .== s,:discharge_efficiency][1]
    ) #TODO: add delta_T
    @constraint(model, SOEDnEvol[s in S, t in T], 
        SOEDN[s,t]  == SOEDN[s,t-1] + CH[s,t]*storage[storage.r_id .== s,:charge_efficiency][1] - (DIS[s,t] + μ_dn*RESUP[s,t])/storage[storage.r_id .== s,:discharge_efficiency][1]
    ) #TODO: add delta_T
    
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

function add_energy_reserve_constraints(model, reserve, loads, gen_df, storage, storage_link_constraint)
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]
    _, G_thermal, _, __, ___, ____ = create_generators_sets(gen_df)
    T, _____ =  create_time_sets(loads)

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
    # (1) Reserves limited by committed capacity of generator
    @constraint(model, EnergyResUpCap[i in G_thermal, j in T, t in T; j <= t],
        ERESUP[i, j, t] <= sum(
            COMMIT[i,tt]*gen_df[gen_df.r_id .== i, :existing_cap_mw][1] - GEN[i,tt] for tt in T if (tt >= j)&(tt <= t)
        )
    )
    @constraint(model, EnergyResDownCap[i in G_thermal, j in T, t in T; j <= t],
        ERESDN[i, j, t] <= sum(
            GEN[i,tt] - COMMIT[i,tt]*gen_df[gen_df.r_id .==i,:existing_cap_mw][1]*gen_df[gen_df.r_id .==i,:min_power][1] for tt in T if (tt >= j)&(tt <= t) #TODO: calculation should be done at input file
        )
    )
    # (2) Reserves limited by ramp rates
    @constraint(model, EnergyResUpRamp[i in G_thermal, j in T, t in T; j <= t],
        ERESUP[i, j, t] <=  sum(
            gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_up_percentage][1] for tt in T if (tt >= j)&(tt <= t) #TODO: calculation should be done at input file
        )
    )
    @constraint(model, EnergyResDnRamp[i in G_thermal, j in T, t in T; j <= t],
        ERESDN[i, j, t] <=  sum(
            gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1] for tt in T if (tt >= j)&(tt <= t) #TODO: calculation should be done at input file
        )
    )

    # @constraint(model, EnergyResUpZero[i in G_thermal, j in T, t in T; j <= t],
    #     ERESUP[i, j, t] == 0
    # )
    # @constraint(model, EnergyResDnZero[i in G_thermal, j in T, t in T; j <= t],
    #     ERESDN[i, j, t] == 0
    # )

    # (3) Storage reserve
    if !isnothing(storage)
        @constraint(model, EnergyResUpStorage[s in S, j in T, t in T; j <= t],
            ERESUP[s, j, t] <= (SOE[s,t]- storage[storage.r_id .== s,:min_energy_mwh][1])*storage[storage.r_id .== s,:discharge_efficiency][1] #TODO: include delta_T
        )
        @constraint(model, EnergyResDownStorage[s in S, j in T, t in T; j <= t],
            ERESDN[s, j, t] <= (storage[storage.r_id .== s,:max_energy_mwh][1] - SOE[s,t])/storage[storage.r_id .== s,:charge_efficiency][1] #TODO: include delta_T
        )
        if storage_link_constraint
            @constraint(model, EnergyResUpLink[s in S, j in T, t in T; j <= t],
            ERESUP[s, j, t] == sum(ERESUP[s, tt, tt] for tt in T if (tt >= j)&(tt <= t))
            )
            @constraint(model, EnergyResDownLink[s in S, j in T, t in T; j <= t],
                ERESDN[s, j, t] == sum(ERESDN[s, tt, tt] for tt in T if (tt >= j)&(tt <= t))
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

function construct_unit_commitment(gen_df, loads, gen_variable; kwargs...)
    storage = nothing
    storage_envelopes = false
    storage_link_constraint =  get(kwargs, :storage_link_constraint, false)
    μ_up = get(kwargs, :μ_up, 1)
    μ_dn = get(kwargs, :μ_dn, 1)
    VRESERVE = get(kwargs, :value_reserve, 1e-6)
    mip_gap = get(kwargs, :mip_gap, 1e-8)

    println("Constructing UC...")
    uc = unit_commitment(gen_df, loads, gen_variable, mip_gap)
    #TODO: modify this part to format arg = get(kwargs, key, default)
    if haskey(kwargs,:storage)
        println("Adding storage...")
        storage = kwargs[:storage]
        if haskey(kwargs,:storage_envelopes)
            if kwargs[:storage_envelopes]
                storage_envelopes = true
            end
        end
        add_storage(uc, storage, loads, gen_df)
    end
    if haskey(kwargs,:ramp_constraints)
        if kwargs[:ramp_constraints] == true
            println("Adding ramp constraints...")   
            add_ramp_constraints(uc, loads, gen_df)
        end
    end
    if haskey(kwargs,:reserve)
        println("Adding reserve constraints...")
        add_reserve_constraints(uc, kwargs[:reserve], loads, gen_df, storage, storage_envelopes, μ_up, μ_dn, VRESERVE)
    end
    if haskey(kwargs,:energy_reserve)
        println("Adding energy reserve constraints...")
        add_energy_reserve_constraints(uc, kwargs[:energy_reserve], loads, gen_df, storage, storage_link_constraint)
    end
    # println("...done")
    return uc
end

function solve_unit_commitment(gen_df, loads, gen_variable; kwargs...)
    uc = construct_unit_commitment(gen_df, loads, gen_variable; kwargs...)
    # relax_integrality(uc)
    # include("./debugging_ignore.jl")
    # set_optimizer_attribute(model, "OutputFlag", 1)
    # @infiltrate
    optimize!(uc)
    solution = get_solution(uc)
    if haskey(kwargs,:enriched_solution)
        if kwargs[:enriched_solution] == true
            return enrich_dfs(solution, gen_df, loads, gen_variable; kwargs...) #TODO: remove kwargs 
        end
    end
    return solution
end
using JuMP
using Gurobi
# include("../utils.jl")
# include("../post_processing.jl")
include("./utils.jl")

function SUC(gen_df, gen_variable, scenarios, mip_gap, VLOL = 10^4, VLGEN = 0)
    # model = direct_model(Gurobi.Optimizer(GRB_ENV ))
    prob = scenarios.probability
    demand = scenarios.demand
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", mip_gap)
    @variable(model, MIPGap in Parameter(mip_gap))
    # set_optimizer_attribute(model, "LogFile", "./output/log_file.txt")
    set_optimizer_attribute(model, "OutputFlag", 1)

    # model = Model(HiGHS.Optimizer)
    # set_optimizer_attribute(model, "mip_rel_gap", mip_gap)
    # G, G_thermal, _, G_var, G_nonvar, G_nt_nonvar = create_generators_sets(gen_df)
    # T, T_red = create_time_sets(demand)
    # Σ = create_scenarios_sets(prob)
    sets = get_sets(gen_df, demand, prob)
    G = sets.G
    G_thermal = sets.G_thermal
    T = sets.T
    Σ = sets.Σ

    VLOL = convert_to_indexed_vector(VLOL, T)
    VLGEN = convert_to_indexed_vector(VLGEN, T)

    @variables(model, begin
        GEN[G,T,Σ] >= 0     # generation
        LOL[T,Σ] >= 0
        LGEN[T,Σ] >= 0
        COMMIT[G_thermal,T,Σ], Bin # commitment status (Bin=binary)
        START[G_thermal,T,Σ], Bin  # startup decision
        SHUT[G_thermal,T,Σ], Bin   # shutdown decision
    end)
              
  # Objective function
      # Sum of variable costs + start-up costs for all generators and time periods
      # TODO: add delta_T
    
    add_OPEX(model, gen_df, sets, scenarios, VLOL, VLGEN)

    @objective(model, Min,  #TODO: move at the end of the constructor
        model[:OPEX]
    )
    # Demand balance constraint (supply must = demand in all time periods)
    @expression(model, SupplyDemand[t in T, σ in Σ],
        sum(GEN[g,t,σ] for g in G) + + LOL[t,σ] - LGEN[t,σ]
    )
    @constraint(model, SupplyDemandBalance[t in T, σ in Σ], 
        SupplyDemand[t,σ] == demand[demand.hour .== t,σ][1]
    )

    add_capacity_constraints(model, gen_df, gen_variable, sets)
    add_commitment_logic(model, gen_df, sets)

    constrain_to_deterministic(model, :GEN, G_thermal)
    constrain_to_deterministic(model, :GEN, sets.G_nt_nonvar)
    constrain_to_deterministic(model, :COMMIT)
    constrain_to_deterministic(model, :START)
    constrain_to_deterministic(model, :SHUT)
    return model
end

function constrain_to_deterministic(model, var_symbol, GG = nothing)
    var = model[var_symbol]
    if isnothing(GG)
        GG = axes(var)[1]
        subset = ""
    else
        subset = "$(minimum(GG))_to_$(maximum(GG))_"
    end
    var_d_symbol = Symbol("$(string(var_symbol))_$(subset)d")
    T = axes(var)[2]
    Σ = axes(var)[3]
    model[var_d_symbol] = @variable(model, [GG,T], base_name = string(var_d_symbol))
    for g in GG, t in T
        # @constraint(model, [σ in Σ], var[g,t,σ] == model[var_d_symbol][g,t]) # for now we declare anonimous constraints
        model[Symbol("$(string(var_d_symbol))_constraint")] = @constraint(model, [σ in Σ],
            var[g,t,σ] == model[var_d_symbol][g,t],
            base_name = "$(string(var_d_symbol))_constraint")
    end
end

function add_OPEX(model, gen_df, sets, scenarios, VLOL, VLGEN)
    prob = scenarios.probability
    G_thermal = sets.G_thermal
    G_nt_nonvar = sets.G_nt_nonvar
    G_nonvar = sets.G_nonvar
    G_var = sets.G_var
    T = sets.T
    Σ = sets.Σ

    GEN = model[:GEN]
    START = model[:START]
    COMMIT = model[:COMMIT]
    LOL = model[:LOL]
    LGEN = model[:LGEN]

    @expression(model, StartupFixedCost[σ in Σ],
        sum(gen_df[gen_df.r_id .== g,:start_cost_per_mw][1]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*START[g,t,σ] for g in G_thermal for t in T) + 
        sum(gen_df[gen_df.r_id .== g,:fixed_om_cost_per_mw_per_hour][1]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*COMMIT[g,t,σ] for g in G_thermal for t in T) + 
        sum(gen_df[gen_df.r_id .== g,:fixed_om_cost_per_mw_per_hour][1]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1] for g in G_nt_nonvar for t in T)
    )

    @expression(model, OperationalCost[σ in Σ],
        sum((gen_df[gen_df.r_id .== g,:heat_rate_mmbtu_per_mwh][1]*gen_df[gen_df.r_id .== g,:fuel_cost][1] + gen_df[gen_df.r_id .== g,:var_om_cost_per_mwh][1])*GEN[g,t,σ] for g in G_nonvar for t in T) +
        sum(gen_df[gen_df.r_id .== g,:var_om_cost_per_mwh][1]*GEN[g,t,σ] for g in G_var for t in T)
    )

    # @expression(model, OperationalCost[σ in Σ],
    #     sum(prob[prob.scenario .== σ,:probability][1]*(model[:ScenarioOperationalCost][σ]) for σ in Σ)
    # )

    @expression(model, LOLCost[σ in Σ],
        sum(VLOL[t]*LOL[t,σ] for t in T)
    )
    # @expression(model, VOLL[σ in Σ],
    #     sum(prob[prob.scenario .== σ,:probability][1]*(model[:ScenarioVOLL][σ]) for σ in Σ)
    # )

    @expression(model, LGENCost[σ in Σ],
        sum(VLGEN[t]*LGEN[t,σ] for t in T)
    )
    # @expression(model, VOLL[σ in Σ],
    #     sum(prob[prob.scenario .== σ,:probability][1]*(model[:ScenarioVLGEN][σ]) for σ in Σ)
    # )

    @expression(model, OPEX,
        sum(prob[prob.scenario .== σ,:probability][1]*(model[:StartupFixedCost][σ] + model[:OperationalCost][σ] + model[:LOLCost][σ] + model[:LGENCost][σ]) for σ in Σ)
    )
    @variable(model, VLOL[t in keys(VLOL)] in Parameter(VLOL[t])) # for post-processing purposes
    @variable(model, VLGEN[t in keys(VLGEN)] in Parameter(VLGEN[t])) # for post-processing purposes
end

function add_capacity_constraints(model, gen_df, gen_variable, sets)
    G_thermal = sets.G_thermal
    G_nt_nonvar = sets.G_nt_nonvar
    T = sets.T
    Σ = sets.Σ
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]
    
    @constraint(model, Cap_thermal_min[g in G_thermal, t in T, σ in Σ], 
        GEN[g,t,σ] >= COMMIT[g,t,σ]*gen_df[gen_df.r_id .== g, :existing_cap_mw][1]*gen_df[gen_df.r_id .== g, :min_power][1] 
    ) 
    @constraint(model, Cap_thermal_max[g in G_thermal, t in T, σ in Σ], 
        GEN[g,t,σ] <= COMMIT[g,t,σ]*gen_df[gen_df.r_id .== g, :existing_cap_mw][1]
    ) 

    # 2. non-variable generation not requiring commitment
    @constraint(model, Cap_nt_nonvar[g in G_nt_nonvar, t in T, σ in Σ], 
        GEN[g,t,σ] <= gen_df[gen_df.r_id .== g, :existing_cap_mw][1]
    )
    # 3. variable generation, accounting for hourly capacity factor
    # TODO: The way this constraint is declared does not follow general style
    # Needs to be redefined at each ED
    @constraint(model, Cap_var[g in 1:nrow(gen_variable), σ in Σ], 
            GEN[gen_variable[g,:r_id], gen_variable[g,:hour], σ] <= 
                        gen_variable[g,:cf] *
                        gen_variable[g,:existing_cap_mw]
                    )
end

function add_commitment_logic(model, gen_df, sets)
    G_thermal = sets.G_thermal
    T = sets.T
    T_red = sets.T_red
    Σ = sets.Σ
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]
    START = model[:START]
    SHUT = model[:SHUT]
     # 1. Minimum up time
     @constraint(model, Startup[g in G_thermal, t in T, σ in Σ],
        COMMIT[g,t,σ] >= sum(START[g,tt,σ] for tt in intersect(T, (t-gen_df[gen_df.r_id .== g,:up_time][1]):t))
    )

    # 2. Minimum down time
    @constraint(model, Shutdown[g in G_thermal, t in T, σ in Σ],
        1-COMMIT[g,t,σ] >= sum(SHUT[g, tt,σ] for tt in intersect(T, (t-gen_df[gen_df.r_id .== g,:down_time][1]):t))
    )

    # 3. Start up/down logic
    @constraint(model, CommitmentStatus[g in G_thermal, t in T_red, σ in Σ],
        COMMIT[g,t+1,σ] - COMMIT[g,t,σ] == START[g,t+1,σ] - SHUT[g,t+1,σ]
    )
end

function add_storage_s(model, storage, scenarios, sets::NamedTuple, expected_min_SOE::Bool)
    prob = scenarios.probability
    demand = scenarios.demand
    
    GEN = model[:GEN]
    T = sets.T
    T_incr = copy(T)
    pushfirst!(T_incr, T_incr[1]-1)
    Σ = sets.Σ
    S = create_storage_sets(storage)

    
    # START = model[:START]
    @variables(model, begin
        CH[S,T,Σ] >= 0
        DIS[S,T,Σ] >= 0
        SOE[S,T_incr,Σ] >= 0 # T_incr captures SOE at t = T[1]-1
        M[S,T,Σ], Bin
    end)

    # Redefinition of objecive function
    @expression(model, StorageOperationalCost[σ in Σ],
        sum(storage[storage.r_id .== s,:var_om_cost_per_mwh][1]*(CH[s,t,σ] + DIS[s,t,σ]) for s in S, t in T)
    )
    # @expression(model, StorageOperationalCost,
    #     sum(model[:ScenarioStorageOperationalCost][σ] for σ in Σ)
    # )

    OPEX = model[:OPEX]
    remove_variable_constraint(model, :OPEX, false)
    @expression(model, OPEX,
        OPEX + sum(prob[prob.scenario .== σ,:probability][1]*model[:StorageOperationalCost][σ] for σ in Σ)
    )
    @objective(model, Min, #TODO: move towars the end
        model[:OPEX]
    )
    
    # Redefinition of supply-demand balance expression and constraint
    SupplyDemand = model[:SupplyDemand]
    unregister(model, :SupplyDemand)
    @expression(model, SupplyDemand[t in T, σ in Σ],
        SupplyDemand[t,σ] - sum(CH[s,t,σ] - DIS[s,t,σ] for s in S)
    )
    SupplyDemandBalance = model[:SupplyDemandBalance]
    delete.(model, SupplyDemandBalance) # Constraints must be deleted also
    unregister(model, :SupplyDemandBalance)
    @constraint(model, SupplyDemandBalance[t in T, σ in Σ], 
        SupplyDemand[t,σ] == demand[demand.hour .== t,σ][1]
    )
    # Charging-discharging logic
    @constraint(model, ChargeLogic[s in S, t in T, σ in Σ],
        CH[s,t,σ] <= storage[storage.r_id .== s,:existing_cap_mw][1]*M[s,t,σ]
    )
    @constraint(model, DischargeLogic[s in S, t in T, σ in Σ],
        DIS[s,t,σ] <= storage[storage.r_id .== s,:existing_cap_mw][1]*(1-M[s,t,σ])
    )
    
    # Storage constraints
    @constraint(model, SOEEvol[s in S, t in T, σ in Σ], 
        SOE[s,t,σ] == SOE[s,t-1,σ] + CH[s,t,σ]*storage[storage.r_id .== s,:charge_efficiency][1] - DIS[s,t,σ]/storage[storage.r_id .== s,:discharge_efficiency][1]
    ) #TODO: add delta_T

    @constraint(model, SOEMax[s in S, t in T, σ in Σ],
        SOE[s,t,σ] <= storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    @constraint(model, SOEMin[s in S, t in T, σ in Σ],
        SOE[s,t,σ] >= storage[storage.r_id .== s,:min_energy_mwh][1]
    )
    @constraint(model, CHMin[s in S, t in T, σ in Σ],
        CH[s,t,σ] >= storage[storage.r_id .== s,:existing_cap_mw][1]*storage[storage.r_id .== s,:min_power][1] #TODO: calculation should be done at input file
    )
    @constraint(model, DISMin[s in S, t in T, σ in Σ],
        DIS[s,t,σ] >= storage[storage.r_id .== s,:existing_cap_mw][1]*storage[storage.r_id .== s,:min_power][1] #TODO: calculation should be done at input file
    )
    # SOE_T_initial = SOE_0
    @constraint(model, SOEO[s in S, σ in Σ], #TODO: replace by T
        SOE[s,T_incr[1],σ] == storage[storage.r_id .== s,:initial_energy_proportion][1]*storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    if expected_min_SOE
        @constraint(model, SOEFinalMin[s in S],
            sum(prob[prob.scenario .== σ,:probability][1]*SOE[s,T[end],σ] for σ in Σ) >= storage[storage.r_id .== s,:final_energy_proportion][1]*storage[storage.r_id .== s,:max_energy_mwh][1]
        )
        # @constraint(model, MaxSOEFinal[s in S],
        #     sum(prob[prob.scenario .== σ,:probability][1]*SOE[s,T[end],σ] for σ in Σ) <= storage[storage.r_id .== s,:initial_energy_proportion][1]*storage[storage.r_id .== s,:max_energy_mwh][1]
        # )
    else
        @constraint(model, SOEFinal[s in S],
        sum(prob[prob.scenario .== σ,:probability][1]*SOE[s,T[end],σ] for σ in Σ) == storage[storage.r_id .== s,:initial_energy_proportion][1]*storage[storage.r_id .== s,:max_energy_mwh][1]
    )
    end
end

function add_ramp_constraints_s(model, gen_df, sets::NamedTuple)
    # G, G_thermal, G_nonthermal, __, ___, ____ = create_generators_sets(gen_df)  
    G_thermal = sets.G_thermal
    G_nonthermal = sets.G_nonthermal
    T= sets.T
    T_red= sets.T_red
    Σ= sets.Σ
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]

    T = axes(GEN)[2]
    T_red = T[1:end-1]
    Σ = axes(GEN)[3]

    # New auxiliary variable GENAUX for generation above the minimum output level
    @variable(model, GENAUX[G_thermal, T, Σ] >= 0)
    
    # for committed thermal units (only created for thermal generators)
    @constraint(model, AuxGen[g in G_thermal, t in T, σ in Σ],
        GENAUX[g,t,σ] == GEN[g,t,σ] - COMMIT[g,t,σ]*gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:min_power][1]
    )
    
    # Ramp equations for thermal generators (constraining GENAUX)
    @constraint(model, RampUp_thermal[g in G_thermal, t in T_red, σ in Σ], 
        GENAUX[g,t+1,σ] - GENAUX[g,t,σ] <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_up_percentage][1]
    )

    @constraint(model, RampDn_thermal[g in G_thermal, t in T_red, σ in Σ], 
        GENAUX[g,t,σ] - GENAUX[g,t+1,σ] <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_dn_percentage][1]
    )

    # Ramp equations for non-thermal generators (constraining total generation GEN)
    @constraint(model, RampUp_nonthermal[g in G_nonthermal, t in T_red, σ in Σ], 
        GEN[g,t+1,σ] - GEN[g,t,σ] <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_up_percentage][1]
    )

    @constraint(model, RampDn_nonthermal[g in G_nonthermal, t in T_red, σ in Σ], 
        GEN[g,t,σ] - GEN[g,t+1,σ] <= gen_df[gen_df.r_id .== g,:existing_cap_mw][1]*gen_df[gen_df.r_id .== g,:ramp_dn_percentage][1]
    )
end









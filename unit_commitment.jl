
using JuMP
using Gurobi
using DataFrames
# using HiGHS
include("./utils.jl")

function value_to_df(var)
    if var isa JuMP.Containers.DenseAxisArray
        return value_to_df_2dim(var)
    elseif var isa JuMP.Containers.SparseAxisArray
        return value_to_df_multidim(var, [:r_id, :hour_i, :hour])
    else
        print("Could not identify type of output. Returnin initial variable...")
        return var
    end
end

function value_to_df_multidim(var, new_index_list = nothing)
    size_of_indices = length(first(keys(value.(var).data)))
    indices = [Symbol("index_"*string(i)) for i in 1:size_of_indices]
    out = DataFrame(index = collect(keys(value.(var).data)), value = collect(values(value.(var).data)))
    transform!(out, :index .=> [ByRow(x -> x[i]) .=> indices[i] for i in 1:size_of_indices])
    select!(out, Not(:index))
    if !isnothing(new_index_list)
        for i in 1:length(new_index_list)
            rename!(out, indices[i] => new_index_list[i])
        end
    end
    return out
end

function value_to_df_2dim(var)
    solution = DataFrame(value.(var).data, :auto)
    ax1 = value.(var).axes[1]
    ax2 = value.(var).axes[2]
    cols = names(solution)
    insertcols!(solution, 1, :r_id => ax1)
    solution = stack(solution, Not(:r_id), variable_name=:hour)
    solution.hour = foldl(replace, [cols[i] => ax2[i] for i in 1:length(ax2)], init=solution.hour)
    #   rename!(solution, :value => :gen)
    solution.hour = convert.(Int64,solution.hour)
    return solution
end

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
    @objective(model, Min, 
        sum( (gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1] * gen_df[gen_df.r_id .== i,:fuel_cost][1] +
            gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1]) * GEN[i,t] 
                        for i in G_nonvar for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1] * GEN[i,t] 
                        for i in G_var for t in T)  + 
        sum(gen_df[gen_df.r_id .== i,:start_cost_per_mw][1] * 
            gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
            START[i,t] 
                        for i in G_thermal for t in T)
    )

    # Demand balance constraint (supply must = demand in all time periods)
    @constraint(model, SupplyDemandBalance[t in T], 
        sum(GEN[i,t] for i in G) == loads[loads.hour .== t,:demand][1])

    # Capacity constraints 
    # 1. thermal generators requiring commitment
    @constraint(model, Cap_thermal_min[i in G_thermal, t in T], 
        GEN[i,t] >= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])
    @constraint(model, Cap_thermal_max[i in G_thermal, t in T], 
        GEN[i,t] <= COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1])

    # 2. non-variable generation not requiring commitment
    @constraint(model, Cap_nt_nonvar[i in G_nt_nonvar, t in T], 
        GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1])

    # 3. variable generation, accounting for hourly capacity factor
    # TODO: The way this constraint is declared does not follow general style
    @constraint(model, Cap_var[i in 1:nrow(gen_variable)], 
            GEN[gen_variable[i,:r_id], gen_variable[i,:hour] ] <= 
                        gen_variable[i,:cf] *
                        gen_variable[i,:existing_cap_mw])

    # Unit commitment constraints
    # 1. Minimum up time
    @constraint(model, Startup[i in G_thermal, t in T],
        COMMIT[i, t] >= sum(START[i, tt] 
                        for tt in intersect(T,
                            (t-gen_df[gen_df.r_id .== i,:up_time][1]):t)))

    # 2. Minimum down time
    @constraint(model, Shutdown[i in G_thermal, t in T],
        1-COMMIT[i, t] >= sum(SHUT[i, tt] 
                        for tt in intersect(T,
                            (t-gen_df[gen_df.r_id .== i,:down_time][1]):t)))

    # 3. Commitment state
    @constraint(model, CommitmentStatus[i in G_thermal, t in T_red],
        COMMIT[i,t+1] - COMMIT[i,t] == START[i,t+1] - SHUT[i,t+1])
    return model
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
        GENAUX[i,t] == GEN[i,t] - COMMIT[i, t] * gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
                        gen_df[gen_df.r_id .== i,:min_power][1])
    
    # Ramp equations for thermal generators (constraining GENAUX)
    @constraint(model, RampUp_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t+1] - GENAUX[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_up_percentage][1])

    @constraint(model, RampDn_thermal[i in G_thermal, t in T_red], 
        GENAUX[i,t] - GENAUX[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1])

    # Ramp equations for non-thermal generators (constraining total generation GEN)
    @constraint(model, RampUp_nonthermal[i in G_nonthermal, t in T_red], 
        GEN[i,t+1] - GEN[i,t] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                                 gen_df[gen_df.r_id .== i,:ramp_up_percentage][1])

    # @constraint(model, RampDn[i in G, t in T_red], 
    #     GEN[i,t] - GEN[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
    #                              gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1])

    @constraint(model, RampDn_nonthermal[i in G_nonthermal, t in T_red], 
    GEN[i,t] - GEN[i,t+1] <= gen_df[gen_df.r_id .== i,:existing_cap_mw][1] * 
                            gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1])
end

function add_reserve_constraints(model, reserve, loads, gen_df)
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]
    _, G_thermal, _, __, ___, ____ = create_generators_sets(gen_df)
    T, _____ =  create_time_sets(loads)
    @variables(model, begin
        RESUP[G_thermal, T] >= 0
        RESDN[G_thermal, T] >= 0
    end)
    # (1) Reserves limited by committed capacity of generator
    @constraint(model, ResUpCap[i in G_thermal, t in T],
        RESUP[i,t] <= COMMIT[i,t]*gen_df[gen_df.r_id .==i,:existing_cap_mw][1] - GEN[i,t])

    @constraint(model, ResDnCap[i in G_thermal, t in T],
        RESDN[i,t] <= GEN[i,t] - COMMIT[i,t]*gen_df[gen_df.r_id .==i,:existing_cap_mw][1]*gen_df[gen_df.r_id .==i,:min_power][1]) #TODO: calculation should be done at input file
    
    # (2) Reserves limited by ramp rates
    @constraint(model, ResUpRamp[i in G_thermal, t in T],
        RESUP[i,t] <=  gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_up_percentage][1] )

    @constraint(model, ResDnRamp[i in G_thermal, t in T],
        RESDN[i,t] <=  gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1] )

    # (3) Overall reserve requirements
    @constraint(model, ResUpRequirement[t in T],
        sum(RESUP[i,t] for i in G_thermal) >= reserve[reserve.hour .== t,:reserve_up_MW][1])
    
    @constraint(model, ResDnRequirement[t in T],
        sum(RESDN[i,t] for i in G_thermal) >= reserve[reserve.hour .== t,:reserve_down_MW][1])
end

function add_energy_reserve_constraints(model, reserve, loads, gen_df)
    GEN = model[:GEN]
    COMMIT = model[:COMMIT]
    _, G_thermal, _, __, ___, ____ = create_generators_sets(gen_df)
    T, _____ =  create_time_sets(loads)
    @variables(model, begin
        RESUP[G_thermal, j in T, t in T; j <= t] >= 0
        RESDN[G_thermal, j in T, t in T; j <= t] >= 0
    end)
    # (1) Reserves limited by committed capacity of generator
    @constraint(model, EnergyResUpCap[i in G_thermal, j in T, t in T; j <= t],
        RESUP[i, j, t] <= sum(
            COMMIT[i,tt]*gen_df[gen_df.r_id .== i, :existing_cap_mw][1] - GEN[i,tt] for tt in T if (tt >= j)&(tt <= t)
        )
    )
    @constraint(model, EnergyResDownCap[i in G_thermal, j in T, t in T; j <= t],
        RESDN[i, j, t] <= sum(
            GEN[i,tt] - COMMIT[i,tt]*gen_df[gen_df.r_id .==i,:existing_cap_mw][1]*gen_df[gen_df.r_id .==i,:min_power][1] for tt in T if (tt >= j)&(tt <= t) #TODO: calculation should be done at input file
        )
    )
    # (2) Reserves limited by ramp rates
    @constraint(model, EnergyResUpRamp[i in G_thermal, j in T, t in T; j <= t],
        RESUP[i, j, t] <=  sum(
            gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_up_percentage][1] for tt in T if (tt >= j)&(tt <= t) #TODO: calculation should be done at input file
        )
    )
    @constraint(model, EnergyResDnRamp[i in G_thermal, j in T, t in T; j <= t],
        RESDN[i, j, t] <=  sum(
            gen_df[gen_df.r_id .== i,:existing_cap_mw][1]*gen_df[gen_df.r_id .== i,:ramp_dn_percentage][1] for tt in T if (tt >= j)&(tt <= t) #TODO: calculation should be done at input file
        )
    )
    # (3) Overall reserve requirements
    @constraint(model, EnergyResUpRequirement[j in T, t in T; j <= t],
        sum(RESUP[i, j, t] for i in G_thermal) >= reserve[(reserve.i_hour .== j).&(reserve.t_hour .== t),:reserve_up_MW][1]
    )
 
    @constraint(model, EnerResDnRequirement[j in T, t in T; j <= t],
        sum(RESDN[i, j, t] for i in G_thermal) >= reserve[(reserve.i_hour .== j).&(reserve.t_hour .== t),:reserve_down_MW][1]
    )

end


function add_storage(model, storage, loads, gen_df)
    G, G_thermal, __, G_var, G_nonvar, ___ = create_generators_sets(gen_df)
    T, ____ =  create_time_sets(loads)
    S = create_storage_sets(storage)
    big_M = 1000
    GEN = model[:GEN]
    START = model[:START]
    @variables(model, begin
        CH[S,T] >= 0
        DIS[S,T] >= 0
        SOE[S,T] >= 0
        M[S,T], Bin
    end)

    # Redefinition of objecive function
    @objective(model, Min, 
        sum((gen_df[gen_df.r_id .== i,:heat_rate_mmbtu_per_mwh][1] * gen_df[gen_df.r_id .== i,:fuel_cost][1] +
            gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1]) * GEN[i,t] 
                        for i in G_nonvar for t in T) + 
        sum(gen_df[gen_df.r_id .== i,:var_om_cost_per_mwh][1] * GEN[i,t] 
                        for i in G_var for t in T)  + 
        sum(gen_df[gen_df.r_id .== i,:start_cost_per_mw][1] * 
            gen_df[gen_df.r_id .== i,:existing_cap_mw][1] *
            START[i,t]
                        for i in G_thermal for t in T) +
        sum(storage[storage.r_id .== s,:var_om_cost_per_mwh][1]*(CH[s,t] + DIS[s,t]) for s in S, t in T) #TODO: add delta_T
    )
    # Redefinition of supply-demand balance constraint
    SupplyDemandBalance = model[:SupplyDemandBalance]
    delete.(model, SupplyDemandBalance)
    unregister(model, :SupplyDemandBalance)
    @constraint(model, SupplyDemandBalance[t in T], 
        sum(GEN[i,t] for i in G) - sum(CH[s,t] - DIS[s,t] for s in S) == loads[loads.hour .== t,:demand][1]
    )
    # Charging-discharging logic
    @constraint(model, ChargeLogic[s in S, t in T],
        CH[s,t] <= big_M*M[s,t]
    )
    @constraint(model, DischargeLogic[s in S, t in T],
        DIS[s,t] <= big_M*(1-M[s,t])
    )

    # Storage constraints
    @constraint(model, SOEEvol[s in S, t in T[2:end-1]], 
        SOE[s,t] - (SOE[s,t-1] + CH[s,t]*storage[storage.r_id .== s,:charge_efficiency][1] - DIS[s,t]/storage[storage.r_id .== s,:discharge_efficiency][1]) == 0
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
    # SOE_T_initial = SOE_T_final
    @constraint(model, SOEO_SOEFinal[s in S, t = [T[1],T[end]]],
        SOE[s,t] == storage[storage.r_id .== s,:initial_energy_proportion][1]*storage[storage.r_id .== s,:max_energy_mwh][1]
    )

end

function get_solution_variables(model)
    variables_to_get = [:GEN, :COMMIT, :SHUT, :CH, :DIS, :SOE, :RESUP, :RESDN]
    return NamedTuple(k => value_to_df(model[k]) for k in intersect(keys(object_dictionary(model)), variables_to_get))
end

function add_reserve(enriched_solution_value, solution)
    # @infiltrate
    reserve = innerjoin(
        rename(solution.RESUP, :value => :reserve_up_MW),
        rename(solution.RESDN, :value => :reserve_down_MW),
        on = [:r_id, :hour]
    )
    # leftjoin allows to add reserve to frame in the left
    return leftjoin(enriched_solution_value, reserve, on = [:r_id, :hour])
end

function to_enriched_df(solution, gen_df, loads, gen_variable; kwargs...)
    #TODO: deal with missing values
    # Curtailment calculation
    curtail = innerjoin(gen_variable, solution.GEN, on = [:r_id, :hour])
    curtail.value = curtail.cf .* curtail.existing_cap_mw - curtail.value
    generation = outerjoin(
        outerjoin(
            rename(solution.GEN, :value => :production_MW),
            rename(curtail[!,[:r_id, :hour, :value]], :value => :curtailment_MW), 
            on = [:r_id, :hour]
        ),
        gen_df[!,[:r_id, :resource, :gen_full]],
        on = :r_id
    )
    out = Dict(:generation => generation)
    if haskey(kwargs, :storage)
        storage = innerjoin(
            rename(solution.CH, :value => :charge_MW),
            rename(solution.DIS, :value => :discharge_MW),
            rename(solution.SOE, :value => :SOE_MWh),
            on = [:r_id, :hour]
        )
        storage = innerjoin(
            storage,
            kwargs[:storage][!,[:r_id, :resource, :storage_full]],
            on = :r_id
        )
        out[:storage] = storage
    end
    if haskey(solution, :RESUP) & haskey(solution, :RESDN)
        for (key,value) in out out[key] = add_reserve(value, solution) end
    end
    demand =  rename(loads, :demand => :demand_MW)
    demand.r_id .= missing
    demand.resource .= "total"
    out[:demand] = demand
    return NamedTuple(out)
end 

function solve_unit_commitment(gen_df, loads, gen_variable, mip_gap; kwargs...) 
    uc = unit_commitment(gen_df, loads, gen_variable, mip_gap)

    if haskey(kwargs,:ramp_constraints)
        if kwargs[:ramp_constraints] == true
            println("Adding ramp constraints...")
            add_ramp_constraints(uc, loads, gen_df)
        end
    end
    if haskey(kwargs,:reserve)
        println("Adding reserve constraints...")
        add_reserve_constraints(uc, kwargs[:reserve], loads, gen_df)
    end
    if haskey(kwargs,:energy_reserve)
        println("Adding energy reserve constraints...")
        add_energy_reserve_constraints(uc, kwargs[:energy_reserve], loads, gen_df)
    end
    if haskey(kwargs,:storage)
        println("Adding storage...")
        add_storage(uc, kwargs[:storage], loads, gen_df)
    end
    optimize!(uc)
    solution = get_solution_variables(uc)
    if haskey(kwargs,:enriched_solution)
        if kwargs[:enriched_solution] == true
            return to_enriched_df(solution, gen_df, loads, gen_variable, ; kwargs...)
        end
    end
    return solution
end
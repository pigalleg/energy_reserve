using JuMP
using Gurobi
using DataFrames
include("./unit_commitment.jl")
include("./post_processing.jl")

COMMIT, START, SHUT, LOL_, RESUP, RESDN, SOEUP, SOEDN, ERESUP, ERESDN, HOUR, GEN, CH, DIS, ResUpRequirement, ResDnRequirement = :COMMIT, :START, :SHUT, :LOL, :RESUP, :RESDN, :SOEUP, :SOEDN, :ERESUP, :ERESDN, :hour, :GEN, :CH, :DIS, :ResUpRequirement, :ResDnRequirement 
ITERATION = :iteration
DEMAND = :demand
NB_ITERATIONS = 10000

# TODO change gen_variable => gen_varialbe_df, loads => loads_df
function construct_economic_dispatch(uc, loads, remove_reserve_constraints, constrain_dispatch, variables_to_constrain, remove_variables_from_objective, VLOL = 1e6, VLGEN = 1e6)
    #TODO: remove loads from arguments
    println("Constructing EC...")
    # Outputs EC by fixing variables of UC
    T, __ = create_time_sets(loads)
    # ed = JuMP.copy(uc)
    # set_optimizer(ed, Gurobi.Optimizer)
    # optimize!(ed)
    ed = uc
    constrain_decision_variables(ed, constrain_dispatch, remove_variables_from_objective, variables_to_constrain)

    # update objective function with LOL term and LGEN
    @variables(ed, begin 
        LOL[T] >= 0
        LGEN[T] >= 0
        end)
    @objective(ed, Min, 
        objective_function(ed) + VLOL*sum(LOL[t] for t in T) + VLGEN*sum(LGEN[t] for t in T)
    )

    # Update supply-demand balance expression
    SupplyDemand = ed[:SupplyDemand]
    remove_variable_constraint(ed, :SupplyDemand, false)
    @expression(ed, SupplyDemand[t in T],
        SupplyDemand[t] + LOL[t] - LGEN[t]
    )
    
    # Remove SOE's circular constraints
    remove_variable_constraint(ed, :SOEFinal)

    # By default, (energy) reserve constraints are removed 
    if remove_reserve_constraints
        remove_energy_and_reserve_constraints(ed)

    end
    println("...done")
    return ed
end

function constrain_decision_variables(model, constrain_dispatch, remove_variables_from_objective, variables_to_constrain = [GEN, CH, DIS], variables_to_fix =  [COMMIT, START, SHUT] )
    function normalize_reserve_variables(res_up_var_value, res_dn_var_value)
        # Normalization to match ResUpRequirement and ResDnRequirement lower bounds
        T = axes(res_up_var_value)[2]
        RESUP_lower_bounds = [constraint_object(model[ResUpRequirement][t]).set.lower for t in T]
        RESDN_lower_bounds = [constraint_object(model[ResDnRequirement][t]).set.lower for t in T]
        res_up_var_value_normalized  = res_up_var_value./sum(Array(res_up_var_value), dims = 1).*transpose(RESUP_lower_bounds)
        res_dn_var_value_normalized  = res_dn_var_value./sum(Array(res_dn_var_value), dims = 1).*transpose(RESDN_lower_bounds)
        return res_up_var_value_normalized, res_dn_var_value_normalized
    end
    variables_to_fix = [(model[var], value.(model[var])) for var in variables_to_fix]
    if constrain_dispatch & haskey(model,RESUP)

        variables_to_constrain =  [(model[var], value.(model[var])) for var in variables_to_constrain]
        res_up_var, res_up_var_value = (model[RESUP], value.(model[RESUP]))
        res_dn_var, res_dn_var_value = (model[RESDN], value.(model[RESDN]))
        # res_up_var_value, res_dn_var_value = normalize_reserve_variables(res_up_var_value, res_dn_var_value)
        constrain_dispatch_variables_according_to_reserve(model, variables_to_constrain,  res_up_var, res_up_var_value, res_dn_var, res_dn_var_value)
    end 
    fix_decision_variables(model, variables_to_fix, remove_variables_from_objective)
end



function constrain_dispatch_variables_according_to_reserve(model, variables_to_constrain,  res_up_var, res_up_var_value, res_dn_var, res_dn_var_value)
    # Dispatch constrained based on the procured reserve at UC stage
    # Function fixes up to three variable types: :GEN, :CH and :DIS
    function get_variable_base_name(variable)
        return Symbol(match(r"([A-z]+)\[", name(first(variable)))[1])
    end

    function constraint_production_variables(var, var_value, res_up_var, res_up_var_value, res_dn_var, res_dn_var_value)
        G = intersect(axes(res_up_var)[1], axes(var)[1])
        T = intersect(axes(res_up_var)[2], axes(var)[2])
        model[Symbol("$(string(get_variable_base_name(var)))$(string(get_variable_base_name(res_up_var)))")] = @constraint(model, [s in G, t in T], 
            var[s,t] <= var_value[s,t] + res_up_var_value[s,t]
        )
        G = intersect(axes(res_dn_var)[1], axes(var)[1])
        T = intersect(axes(res_dn_var)[2], axes(var)[2])
        model[Symbol("$(string(get_variable_base_name(var)))$(string(get_variable_base_name(res_dn_var)))")] = @constraint(model, [s in G, t in T], 
            -var[s,t] <= -var_value[s,t] + res_dn_var_value[s,t]
        )
        # By default, units not offering reserve will have their dispatch fixed.
        G_to_fix = setdiff(axes(var)[1], G)
        for key in collect(keys(var)) if key.I[1] in G_to_fix
                fix(var[key], var_value[key], force = true) # force is needed becase the variable has bounds defined.
            end
        end
    end

    println("Constraining dispatch to procured reserve...")
    # res_up_var, res_up_var_value = res_up_variables
    # res_dn_var, res_dn_var_value = res_dn_variables
    for (var, var_value) in variables_to_constrain
        # base_name = name(first(var))
        # base_name = Symbol(match(r"([A-z]+)\[", base_name)[1])
        if get_variable_base_name(var) in [GEN, DIS]
            constraint_production_variables(var, var_value, res_up_var, res_up_var_value, res_dn_var, res_dn_var_value)
        elseif get_variable_base_name(var) in [CH]
            # Inversing arugement allows to constraint consumption variables
            constraint_production_variables(var, var_value, res_dn_var, res_dn_var_value, res_up_var, res_up_var_value)
        end
    end
end

function fix_decision_variables(model, variables, remove_variables_from_objective = true)
    println("Fixing decision variables...")
    for (var, var_value) in variables
    # for (var, var_value) in [(model[var], value.(model[var])) for var in decision_variables]
        for key in collect(keys(var))
           fix(var[key], var_value[key]; force = !is_binary(var[key]))
           if remove_variables_from_objective 
                println("Removing fixed decision variables from objective...")
                set_objective_coefficient(model, var[key], 0) 
            end # when set to zero they are removed from objective function
        end
    end
end

function remove_energy_and_reserve_constraints(model)
    println("Removing reserve, energy reserve and envelope constraints...")
    # Remove reserve, energy reserve and storge envelope's associated variables/constraints
    keys = [:ResUpCap, :ResDnCap, :ResUpRamp, :ResDnRamp, :ResUpStorage, :ResDownStorage, :ResUpRequirement, :ResDnRequirement, :SOEUpEvol, :SOEDnEvol, :SOEUP_0, :SOEDN_0, :SOEUPMax, :SOEDNMax, :SOEUPMin, :SOEDNMin, :EnergyResUpCap, :EnergyResDownCap, :EnergyResUpRamp, :EnergyResDnRamp, :EnrgyResUpStorage, :EnergyResDownStorage, :EnergyResUpRequirement, :EnerResDnRequirement]
    for k in keys
        if haskey(model, k)
            remove_variable_constraint(model, k)
        end
    end
end

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

function update_demand(model, loads, key = DEMAND)
    # Update demand values and introduces LOL at supply-demand balance
    T, __ = create_time_sets(loads)
    LOL = model[LOL_]

    if haskey(model, :LOLMax) remove_variable_constraint(model, :LOLMax) end
    @constraint(model, LOLMax[t in T],
        LOL[t]<= loads[loads.hour .== t, key][1]    
    )

    SupplyDemand = model[:SupplyDemand]
    remove_variable_constraint(model, :SupplyDemandBalance)
    @constraint(model, SupplyDemandBalance[t in T], 
        SupplyDemand[t] == loads[loads.hour .== t, key][1]
    )
end

function update_generation(model, gen_variable)
    remove_variable_constraint(model, :Cap_var)
    GEN = model[:GEN]
    @constraint(model, Cap_var[i in 1:nrow(gen_variable)], 
            GEN[gen_variable[i,:r_id], gen_variable[i,:hour] ] <= gen_variable[i,:cf]*gen_variable[i,:existing_cap_mw]
    )
end



function merge_solutions(solutions::Dict, merge_keys = [ITERATION])
    #TODO can be done more elegantly
    # solution_keys = keys(solutions[collect(keys(solutions))[1]]) # assumes all solutions have the same set of keys
    solution_keys = union([keys(v) for (k,v) in solutions]...)
    aux = Dict(k => [] for k in solution_keys)

    for d in keys(solutions), k in intersect(keys(solutions[d]), solution_keys)
        aux_ = DataFrame(collect(repeat([isa(d,Tuple) ? d : tuple(d)], size(solutions[d][k],1))), merge_keys)
        push!(aux[k], hcat(solutions[d][k], aux_))
    end
    return NamedTuple(k => vcat(aux[k]..., cols = :union) for k in keys(aux))
end


function solve_economic_dispatch_(ed, gen_df, loads, gen_variable; kwargs...)
    print("Solving ED...")
    optimize!(ed)
    solution = get_solution(ed)
    if haskey(kwargs,:enriched_solution)
        if kwargs[:enriched_solution] == true
            println("done")
            return enrich_dfs(solution, gen_df, loads, gen_variable; kwargs...) #TODO: remove kwargs 
        end
    end
    println("done")
    return solution
end

function solve_economic_dispatch(gen_df, loads, gen_variable, mip_gap; kwargs...)
    # Parsing arguments...
    remove_reserve_constraints =get(kwargs, :remove_reserve_constraints, true)
    max_iterations = (get(kwargs, :max_iterations, NB_ITERATIONS))
    constrain_dispatch = get(kwargs, :constrain_dispatch, true)
    variables_to_constrain = get(kwargs, :variables_to_constrain, [GEN, CH, DIS])
    remove_variables_from_objective = get(kwargs, :remove_variables_from_objective, false)
    VLOL = get(kwargs, :VLOL, 1e6)
    VLGEN = get(kwargs, :VLGEN, 1e6)
    # parsing end

    uc = construct_unit_commitment(gen_df, loads[!,[HOUR, DEMAND]], gen_variable, mip_gap; kwargs...)
    optimize!(uc)
    ed = construct_economic_dispatch(uc, loads[!,[HOUR, DEMAND]], remove_reserve_constraints, constrain_dispatch, variables_to_constrain, remove_variables_from_objective, VLOL, VLGEN)

    solutions = Dict()
    # for k in collect(propertynames(loads[!, Not(HOUR)]))
    for k in first(propertynames(loads[!, Not(HOUR)]), max_iterations)
        println("")
        println("Montecarlo iteration: $(k)")
        gen_df_k, loads_df_k, gen_variable_k = pre_process_load_gen_variable(gen_df, rename(loads[!,[HOUR,k]], k=>DEMAND), gen_variable)
        update_demand(ed, loads_df_k)
        update_generation(ed, gen_variable_k)
        solutions[k] = solve_economic_dispatch_(ed, gen_df_k, loads_df_k, gen_variable_k; kwargs...)
    end
    return merge_solutions(solutions)
end
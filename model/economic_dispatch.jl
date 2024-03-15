using JuMP
using Gurobi
using DataFrames
include("./unit_commitment.jl")
include("./post_processing.jl")

COMMIT, START, SHUT, LOL_, RESUP, RESDN, SOEUP, SOEDN, ERESUP, ERESDN, HOUR = :COMMIT, :START, :SHUT, :LOL, :RESUP, :RESDN, :SOEUP, :SOEDN, :ERESUP, :ERESDN, :hour 
ITERATION = :iteration
DEMAND = :demand
NB_ITERATIONS = 10000

function construct_economic_dispatch(uc, loads, remove_reserve_constraints = true, VLOL = 10^6)
    #TODO: remove loads from arguments
    println("Constructing EC...")
    # Outputs EC by fixing variables of UC
    T, __ = create_time_sets(loads)
    # ec = copy_model(uc) # TODO: how to copy model?
    ed = uc
    fix_decision_variables(ed)

    # update objective function with LOL term
    @variables(ed, begin LOL[T] >= 0 end)
    @objective(ed, Min, 
        objective_function(ed) + VLOL*sum(LOL[t] for t in T)
    )
    # By default, (energy) reserve constraints are removed 
    if remove_reserve_constraints
        remove_energy_and_reserve_constraints(ed)
    end
    return ed
end

function remove_energy_and_reserve_constraints(model)
    println("Removing reserve, energy reserve and envelope constraints...")
    # Remove reserve, energy reserve and storge envelope's associated variables/constraints
    # keys = [RESUP, RESDN, :ResUpCap, :ResDnCap, :ResUpRamp, :ResDnRamp, :ResUpStorage, :ResDownStorage, :ResUpRequirement, :ResDnRequirement, :SOEUP, :SOEDN, :SOEUpEvol, :SOEDnEvol, :SOEUP_0, :SOEDN_0, :SOEUPMax, :SOEDNMax, :SOEUPMin, :SOEDNMin, ERESUP, ERESDN, :EnergyResUpCap, :EnergyResDownCap, :EnergyResUpRamp, :EnergyResDnRamp, :EnrgyResUpStorage, :EnergyResDownStorage, :EnergyResUpRequirement, :EnerResDnRequirement]
    keys = [:ResUpCap, :ResDnCap, :ResUpRamp, :ResDnRamp, :ResUpStorage, :ResDownStorage, :ResUpRequirement, :ResDnRequirement, :SOEUpEvol, :SOEDnEvol, :SOEUP_0, :SOEDN_0, :SOEUPMax, :SOEDNMax, :SOEUPMin, :SOEDNMin, :EnergyResUpCap, :EnergyResDownCap, :EnergyResUpRamp, :EnergyResDnRamp, :EnrgyResUpStorage, :EnergyResDownStorage, :EnergyResUpRequirement, :EnerResDnRequirement]

    for k in keys
        if haskey(model, k)
            remove_variable_constraint(model, k)
        end
    end
end

function remove_variable_constraint(model, key, delete_ = true)
    # Applies for constraints and variables
    println(key)
    if delete_ delete.(model, model[key]) end # Constraints must be deleted also
    unregister(model, key)
end

function update_demand(model, loads, key)
    # Update demand values and introduces LOL at supply-demand balance
    T, __ = create_time_sets(loads)
    LOL = model[LOL_]

    if haskey(model, :LOLMax) remove_variable_constraint(model, :LOLMax) end
    @constraint(model, LOLMax[t in T],
        LOL[t]<= loads[loads.hour .== t, key][1]    
    )

    SupplyDemand = model[:SupplyDemand]
    remove_variable_constraint(model, :SupplyDemand, false)
    @expression(model, SupplyDemand[t in T],
        SupplyDemand[t] + LOL[t]
    )

    remove_variable_constraint(model, :SupplyDemandBalance)
    @constraint(model, SupplyDemandBalance[t in T], 
        SupplyDemand[t] == loads[loads.hour .== t, key][1]
    )
end

function fix_decision_variables(model, decision_variables::Array = [COMMIT, START, SHUT])
    println("Fixing decision variables...")
    for (var, var_value) in [(model[var], value.(model[var])) for var in decision_variables]
        for key in collect(keys(var))
           fix(var[key], var_value[key]; force = !is_binary(var[key]))
        end
    end
end

function merge_solutions(solutions::Dict, merge_key = ITERATION)
    #TODO can be done more elegantly
    # solution_keys = keys(solutions[collect(keys(solutions))[1]]) # assumes all solutions have the same set of keys
    solution_keys = union([keys(v) for (k,v) in solutions]...)
    out = Dict(k => [] for k in solution_keys)
    for d in keys(solutions), k in intersect(keys(solutions[d]), solution_keys)
        solutions[d][k][!, merge_key] .= d
        push!(out[k], solutions[d][k])
    end
    return NamedTuple(k => vcat(out[k]..., cols = :union) for k in keys(out))
end


function solve_economic_dispatch_(ed, gen_df, loads, gen_variable; kwargs...)
    println("Solving EC...")
    optimize!(ed)
    solution = get_solution(ed)
    if haskey(kwargs,:enriched_solution)
        if kwargs[:enriched_solution] == true
            return enrich_dfs(solution, gen_df, loads, gen_variable; kwargs...) #TODO: remove kwargs 
        end
    end
    return solution
end

function solve_economic_dispatch(gen_df, loads, gen_variable, mip_gap; kwargs...)
    # Parsing arguments...
    remove_reserve_constraints = (get(kwargs, :remove_reserve_constraints, nothing) == true)
    max_iterations = (get(kwargs, :max_iterations, NB_ITERATIONS))
    # parsing end

    uc = construct_unit_commitment(gen_df, loads[!,[HOUR, DEMAND]], gen_variable, mip_gap; kwargs...)
    optimize!(uc)
    ed  = construct_economic_dispatch(uc, loads[!,[HOUR, DEMAND]], remove_reserve_constraints)

    solutions = Dict()
    # for k in collect(propertynames(loads[!, Not(HOUR)]))
    for k in first(propertynames(loads[!, Not(HOUR)]), max_iterations)
        update_demand(ed, loads[!,[HOUR, k]], k)
        solutions[k] = solve_economic_dispatch_(ed, gen_df, loads[!,[HOUR,k]], gen_variable; kwargs...)
    end
    return merge_solutions(solutions)
end

# function solve_economic_dispatch_multiple_demand(multiple_loads, gen_df, loads, gen_variable, mip_gap; kwargs...)
#     uc = construct_unit_commitment(gen_df, loads, gen_variable, mip_gap; kwargs...)
#     optimize!(uc)
#     remove_reserve_constraints = false
#     if haskey(kwargs, :remove_reserve_constraints)
#         remove_reserve_constraints = kwargs[:remove_reserve_constraints] == true
#     end
#     ed  = construct_economic_dispatch(uc, loads, remove_reserve_constraints)
#     @infiltrate
#     update_demand(ed, loads)
#     return solve_economic_dispatch(ed, gen_df, loads, gen_variable; kwargs...)
# end
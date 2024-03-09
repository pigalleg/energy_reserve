using JuMP
using Gurobi
using DataFrames
include("./unit_commitment.jl")
include("./post_processing.jl")

COMMIT, START, SHUT, LOL_ = :COMMIT, :START, :SHUT, :LOL
VLOL = 10^6

function construct_economic_dispatch(uc, loads, remove_reserve_constraints = true)
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

    # Update demand
    update_demand(ed, loads)    

    # By default, (energy) reserve constraints are removed 
    if remove_reserve_constraints
        nothing
    end
    return ed
end

function update_demand(model, loads)
    # Update demand values and introduces LOL on balance
    T, __ = create_time_sets(loads)
    LOL = model[LOL_]
    @constraint(model, LOLMax[t in T],
        LOL[t]<= loads[loads.hour .== t,:demand][1]    
    )

    SupplyDemand = model[:SupplyDemand]
    unregister(model, :SupplyDemand)
    @expression(model, SupplyDemand[t in T],
        SupplyDemand[t] + LOL[t]
    )

    SupplyDemandBalance = model[:SupplyDemandBalance]
    delete.(model, SupplyDemandBalance) # Constraints must be deleted also
    unregister(model, :SupplyDemandBalance)
    @constraint(model, SupplyDemandBalance[t in T], 
        SupplyDemand[t] == 0
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


function solve_economic_dispatch(ed, gen_df, loads, gen_variable; kwargs...)
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

function solve_economic_dispatch_single_demand(gen_df, loads, gen_variable, mip_gap; kwargs...)
    println("Constructing UC...")
    uc = construct_unit_commitment(gen_df, loads, gen_variable, mip_gap; kwargs...)
    optimize!(uc)
    println("Constructing EC...")
    ed  = construct_economic_dispatch(uc, loads)
    return solve_economic_dispatch(ed, gen_df, loads, gen_variable; kwargs...)
end
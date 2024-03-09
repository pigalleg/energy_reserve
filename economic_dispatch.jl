using JuMP
using Gurobi
using DataFrames
include("./unit_commitment.jl")
include("./post_processing.jl")

COMMIT, START, SHUT = :COMMIT, :START, :SHUT
VLOL = 10^6

function construct_economic_dispatch(uc, loads)
    # Outputs EC by fixing variables of UC
    # ec = copy_model(uc)
    # TODO: update objective function here!
    ed = uc
    fix_decision_variables(ed)

    # update objective function with LOL
    
    # remove reserve constraints

    # remove energy constraints
    return ed
end

function update_demand(model, loads)
    # Update demand and introduces LOL
    # TODO: rewrite balance equation every time demand is updated. Is this really needed?
    # @variables(model, begin 
    #     LOL[T] >=0 
    # end)
    # GEN = model[:GEN]
    # u = model[:u]
    # balance = model[:balance]
    # delete.(model, balance)
    # unregister(model, :balance)
    # @constraint(model, balance[t in T], sum(u[s,t]*GEN[s,t] for s in S) - CUR[t]- (D[t]-LOL[t])== 0)
    # @objective(model, Min, sum(sum(u[s,t]*GEN[s,t]*C[s] for s in S) + LOL[t]*VLOLL + CUR[t]*VCUR  for t in T))
    nothing
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
using JuMP
using Gurobi
using DataFrames
include("./unit_commitment.jl")
include("./post_processing.jl")

COMMIT, START, SHUT, LOL_, RESUP, RESDN, SOEUP, SOEDN, ERESUP, ERESDN = :COMMIT, :START, :SHUT, :LOL, :RESUP, :RESDN, :SOEUP, :SOEDN, :ERESUP, :ERESDN

function construct_economic_dispatch(uc, loads, remove_reserve_constraints = true, VLOL = 10^6)
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

    # removed = []
    for k in keys
        if haskey(model, k)
            # push!(removed, k)
            remove_variable_constraint(model, k)
        end
    end
    # @infiltrate
    # display(removed)
    # exit()
end

function remove_variable_constraint(model, key)
    # Applies for constraints and variables
    println(key)
    delete.(model, model[key]) # Constraints must be deleted also
    unregister(model, key)

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
    remove_reserve_constraints = false
    if haskey(kwargs, :remove_reserve_constraints)
        remove_reserve_constraints = kwargs[:remove_reserve_constraints] == true
    end
    ed  = construct_economic_dispatch(uc, loads, remove_reserve_constraints)
    update_demand(ed, loads)
    return solve_economic_dispatch(ed, gen_df, loads, gen_variable; kwargs...)
end
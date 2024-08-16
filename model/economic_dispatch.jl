using JuMP
using Gurobi
using DataFrames
# using Revise
include("./unit_commitment/unit_commitment.jl")
include("./post_processing.jl")
include("../debugging_ignore.jl")
# __revise_mode__ = :eval
COMMIT, START, SHUT, LOL_, RESUP, RESDN, SOEUP, SOEDN, ERESUP, ERESDN, HOUR, GEN, CH, DIS, ResUpRequirement, ResDnRequirement = :COMMIT, :START, :SHUT, :LOL, :RESUP, :RESDN, :SOEUP, :SOEDN, :ERESUP, :ERESDN, :hour, :GEN, :CH, :DIS, :ResUpRequirement, :ResDnRequirement 
ITERATION = :iteration
DEMAND = :demand
NB_ITERATIONS = 10000

function get_multipliers(model)
    CH = model[:CH]
    DIS = model[:DIS]
    RESDNCH = model[:RESDNCH]
    RESUPCH = model[:RESUPCH]
    S = axes(CH)[1]
    T = axes(CH)[2]

    SOE_constraint_list = [constraint_object.(model[:SOEEvol][s,T[1]]).func for s in S]
    η_ch =-map(coefficient, SOE_constraint_list, Array(CH[:,T[1]]))
    
    # We assume that μ=μ(t), but independent of storage unit. We use therefore the first storage to determine the value: η_ch[1]
    SOEUp_constraint_list  = [constraint_object.(model[:SOEUpEvol][S[1],t]).func for t in T]
    μ_dn = -map(coefficient, SOEUp_constraint_list, Array(RESDNCH[S[1],:]))/η_ch[1]
    SOEDN_constraint_list  = [constraint_object.(model[:SOEDnEvol][S[1],t]).func for t in T]
    μ_up =map(coefficient, SOEDN_constraint_list, Array(RESUPCH[S[1],:]))./η_ch[1]
    
    return μ_up, μ_dn
end

# TODO change gen_variable => gen_varialbe_df, loads => loads_df
function construct_economic_dispatch(uc, loads, constraint_dispatch_by_multipliers::Bool, constrain_SOE_by_envelopes::Bool, remove_reserve_constraints::Bool, constrain_dispatch::Bool, constrain_dispatch_by_energy::Bool, bidirectional_storage_reserve::Bool, remove_variables_from_objective::Bool, variables_to_constrain::Vector{Symbol}, VLOL::Union{Float64,Int64,Vector}, VLGEN::Union{Float64,Int64,Vector})
    #TODO: remove loads from arguments
    println("Constructing EC...")
    # Outputs EC by fixing variables of UC
    T, __ = create_time_sets(loads)
    VLOL = convert_to_indexed_vector(VLOL, T)
    VLGEN = convert_to_indexed_vector(VLGEN, T)
    
    # ed, reference_map = copy_model(uc)
    # ed = JuMP.copy(uc)
    # set_optimizer(ed, Gurobi.Optimizer)
    # set_optimizer_attribute(ed, "OutputFlag", 0)
    # set_optimizer_attribute(ed, "MIPGap", get_optimizer_attribute(uc,"MIPGap"))
    # optimize!(ed)
    ed = uc # pointer, uc object will change
    constrain_decision_variables(ed, constraint_dispatch_by_multipliers, constrain_SOE_by_envelopes, constrain_dispatch, constrain_dispatch_by_energy, bidirectional_storage_reserve, remove_variables_from_objective, variables_to_constrain)
    constraint_SOE_final_to_envelopes(ed) # redundant when constrain_SOE_by_envelopes == true
    # update objective function with LOL term and LGEN
    @variables(ed, begin 
        LOL[T] >= 0
        LGEN[T] >= 0
        end)
    @objective(ed, Min, 
        objective_function(ed) + sum(LOL[t]*VLOL[t] + LGEN[t]*VLGEN[t] for t in T)
    )
    # @constraint(ed,
    #     sum(LOL[t] for t in T) == 0 
    # )
    # Update supply-demand balance expression
    # Update of SupplyDemand constraint is performed within the solve_economic_dispatch's loop
    SupplyDemand = ed[:SupplyDemand]
    remove_variable_constraint(ed, :SupplyDemand, false)
    @expression(ed, SupplyDemand[t in T],
        SupplyDemand[t] + LOL[t] - LGEN[t]
    )

    # Remove SOE's circular constraints
    # remove_variable_constraint(ed, :SOEFinal)

    # By default, (energy) reserve constraints are removed 
    if remove_reserve_constraints
        remove_energy_and_reserve_constraints(ed)
    end
    println("...done")
    return ed
end


function constrain_decision_variables(model, constraint_dispatch_by_multipliers::Bool, constrain_SOE_by_envelopes::Bool, constrain_dispatch::Bool, constrain_dispatch_by_energy::Bool, bidirectional_storage_reserve::Bool, remove_variables_from_objective::Bool, variables_to_constrain = [GEN, CH, DIS], variables_to_fix =  [COMMIT, START, SHUT,:RESUP, :RESDN])
    #TODO: split this in multiple functions. Too many arguments.
    variables_to_fix = [(model[var], value.(model[var])) for var in variables_to_fix]
    SOEUP_value, SOEDN_value = generate_envelopes(model)
    if constrain_dispatch & haskey(model,RESUP)
        variables_to_constrain =  [(model[var], value.(model[var])) for var in variables_to_constrain]
        μ_up, μ_dn = get_multipliers(model)
        if !constraint_dispatch_by_multipliers
            μ_up = ones(1, length(μ_up))
            μ_dn = ones(1, length(μ_dn))
        end
        kwargs = Dict(
            :res_up_var => model[:RESUP],
            :res_up_var_value => value.(model[:RESUP]),
            :res_dn_var => model[:RESDN],
            :res_dn_var_value => value.(model[:RESDN]),
            :res_up_ch_var =>  model[:RESUPCH],
            :res_up_ch_var_value => value.(model[:RESUPCH]).*μ_up',
            :res_up_dis_var =>  model[:RESUPDIS],
            :res_up_dis_var_value => value.(model[:RESUPDIS]).*μ_up',
            :res_dn_ch_var =>  model[:RESDNCH],
            :res_dn_ch_var_value => value.(model[:RESDNCH]).*μ_dn',
            :res_dn_dis_var =>  model[:RESDNDIS],
            :res_dn_dis_var_value => value.(model[:RESDNDIS]).*μ_dn',
            )
        # variable extractions must be executeed before modification of the model
        constrain_dispatch_variables_according_to_reserve(model, constrain_dispatch_by_energy, bidirectional_storage_reserve, variables_to_constrain; kwargs...)
    end 
    if constrain_SOE_by_envelopes
        constrain_SOE_to_envelopes(model, SOEUP_value, SOEDN_value)
    end
    @expression(model, SOEUP_EC, SOEUP_value) # added to recover them as output 
    @expression(model, SOEDN_EC, SOEDN_value)
    fix_decision_variables(model, variables_to_fix, remove_variables_from_objective)
end

function constrain_dispatch_variables_according_to_reserve(model, constrain_dispatch_by_energy, bidirectional_storage_reserve, variables_to_constrain; kwargs...)
    # Dispatch constrained based on the procured reserve at UC stage
    # Function fixes up to three variable types: :GEN, :CH and :DIS
    function get_variable_base_name(variable)
        return Symbol(match(r"([A-z]+)\[", name(first(variable)))[1])
    end

    function constraint_production_variables(var, var_value, res_up_var, res_up_var_value; constrain_by_energy = false, lower_bound = false)
        G = intersect(axes(res_up_var)[1], axes(var)[1])
        T = intersect(axes(res_up_var)[2], axes(var)[2])
        name = Symbol("$(string(get_variable_base_name(var)))$(string(get_variable_base_name(res_up_var)))")
        c = !lower_bound ? 1 : -1
        # model[name] = @constraint(model, [s in G, t in T], 
        #     c*var[s,t] <= c*var_value[s,t] + res_up_var_value[s,t]
        #     )
        if !constrain_by_energy
            model[name] = @constraint(model, [s in G, t in T], 
                c*var[s,t] <= c*var_value[s,t] + res_up_var_value[s,t]
            )
            for s in G, t in T # important to identify constraints for debugging
                set_name(model[name][s,t], string(name)*"[$s,$t]")
            end
        else #deprecated
            model[name] = @constraint(model, [s in G], 
                sum(c*var[s,t] for t in T) <= sum(c*var_value[s,t] + res_up_var_value[s,t] for t in T)
            )
            for s in G # important to identify constraints for debugging
                set_name(model[name][s], string(name)*"[$s]")
            end
        end
        # By default, units not offering reserve will have their dispatch fixed.
        G_to_fix = setdiff(axes(var)[1], G)
        for key in collect(keys(var)) if key.I[1] in G_to_fix
                fix(var[key], var_value[key], force = true) # force is needed because the variable has bounds defined.
            end
        end
    end
    println("Constraining dispatch to procured reserve...")
    if bidirectional_storage_reserve
        gen_logic_group = [GEN]
        dis_logic_group = [DIS]
        ch_logic_group = [CH]
        ch_logic_group_2 = []
    else
        gen_logic_group = [GEN, DIS]
        dis_logic_group = []
        ch_logic_group = []
        ch_logic_group_2 = [CH]
    end
    for (var, var_value) in variables_to_constrain
        if get_variable_base_name(var) in gen_logic_group
            # constraint_production_variables(var, var_value, kwargs[:res_up_var], kwargs[:res_up_var_value], kwargs[:res_dn_var], kwargs[:res_dn_var_value])
            # name_up = Symbol("$(string(get_variable_base_name(var)))$(string(get_variable_base_name(res_up_var)))")
            constraint_production_variables(var, var_value, kwargs[:res_up_var], kwargs[:res_up_var_value], constrain_by_energy = constrain_dispatch_by_energy)
            constraint_production_variables(var, var_value, kwargs[:res_dn_var], kwargs[:res_dn_var_value], constrain_by_energy = constrain_dispatch_by_energy, lower_bound = true)
        elseif get_variable_base_name(var) in dis_logic_group
            # constraint_production_variables(var, var_value, kwargs[:res_up_dis_var], kwargs[:res_up_dis_var_value], kwargs[:res_dn_dis_var], kwargs[:res_dn_dis_var_value])
            constraint_production_variables(var, var_value, kwargs[:res_up_dis_var], kwargs[:res_up_dis_var_value], constrain_by_energy = constrain_dispatch_by_energy)
            constraint_production_variables(var, var_value, kwargs[:res_dn_dis_var], kwargs[:res_dn_dis_var_value], constrain_by_energy = constrain_dispatch_by_energy, lower_bound = true)
        elseif get_variable_base_name(var) in ch_logic_group
            # constraint_production_variables(var, var_value, kwargs[:res_dn_ch_var], kwargs[:res_dn_ch_var_value], kwargs[:res_up_ch_var], kwargs[:res_up_ch_var_value])
            constraint_production_variables(var, var_value, kwargs[:res_dn_ch_var], kwargs[:res_dn_ch_var_value], constrain_by_energy = constrain_dispatch_by_energy)
            constraint_production_variables(var, var_value, kwargs[:res_up_ch_var], kwargs[:res_up_ch_var_value], constrain_by_energy = constrain_dispatch_by_energy, lower_bound = true)

        elseif get_variable_base_name(var) in ch_logic_group_2
            # constraint_production_variables(var, var_value, kwargs[:res_dn_var], kwargs[:res_dn_var_value], kwargs[:res_up_var], kwargs[:res_up_var_value])
            constraint_production_variables(var, var_value, kwargs[:res_dn_var], kwargs[:res_dn_var_value], constrain_by_energy = constrain_dispatch_by_energy)
            constraint_production_variables(var, var_value, kwargs[:res_up_var], kwargs[:res_up_var_value], constrain_by_energy = constrain_dispatch_by_energy, lower_bound = true)
        end
    end
end


function generate_envelopes(model)
    # SOEUP_value = value.(model[:SOEUP])
    # SOEDN_value = value.(model[:SOEDN])
    CH = model[:CH]
    DIS = model[:DIS]
    RESDNCH = model[:RESDNCH]
    RESDNDIS = model[:RESDNDIS]
    RESUPCH = model[:RESUPCH]
    RESUPDIS = model[:RESUPDIS]
    S = axes(CH)[1]
    T = axes(CH)[2]
    SOE = model[:SOE]
    T_incr = axes(SOE)[2]

    μ_up, μ_dn = get_multipliers(model)
    SOE_constraint_list = [constraint_object.(model[:SOEEvol][s,T[1]]).func for s in S]
    η_ch =-map(coefficient, SOE_constraint_list, Array(CH[:,T[1]]))
    inv_η_dis = map(coefficient, SOE_constraint_list, Array(DIS[:,T[1]]))

    # SOEPUP_value = +Array(value.(CH)).*η_ch + (Array(value.(RESDNCH)).*η_ch + Array(value.(RESDNDIS)).*inv_η_dis).* μ_dn' # approach 3
    # SOEPUP_value = hcat(zeros(1,size(SOEPUP_value)[1])', SOEPUP_value) # For T[1]-1 no reserves are activated
    # SOEPUP_value = [value(SOE[s,T_incr[1]]) for s in S, t in T_incr] + cumsum(SOEPUP_value; dims = 2)
    
    SOEPUP_value = (Array(value.(RESDNCH)).*η_ch + Array(value.(RESDNDIS)).*inv_η_dis).* μ_dn' #approach 1&2
    SOEPUP_value = hcat(zeros(1,size(SOEPUP_value)[1])', SOEPUP_value) #  #approach 1&2
    SOEPUP_value = Array(value.(SOE))  + cumsum(SOEPUP_value; dims = 2) # approach 1
    # SOEPUP_value = [value(SOE[s,T_incr[1]]) for s in S, t in T_incr] + cumsum(SOEPUP_value; dims = 2) # approach 2


    # SOEPDN_value = -Array(value.(DIS)).*inv_η_dis -(Array(value.(RESUPCH)).*η_ch + Array(value.(RESUPDIS)).*inv_η_dis).* μ_up'# approach 3
    # SOEPDN_value = hcat(zeros(1,size(SOEPDN_value)[1])', SOEPDN_value)
    # SOEPDN_value = [value(SOE[s,T_incr[1]]) for s in S, t in T_incr] + cumsum(SOEPDN_value; dims = 2)

    SOEPDN_value = -(Array(value.(RESUPCH)).*η_ch + Array(value.(RESUPDIS)).*inv_η_dis).* μ_up' #approach 1&2
    SOEPDN_value = hcat(zeros(1,size(SOEPDN_value)[1])', SOEPDN_value) # approach 1&2
    SOEPDN_value = Array(value.(SOE))  + cumsum(SOEPDN_value; dims = 2) # approach 1
    # SOEPDN_value = [value(SOE[s,T_incr[1]]) for s in S, t in T_incr] + cumsum(SOEPDN_value; dims = 2) # approach 2

    SOEMax_value = hcat(Array(normalized_rhs.(model[:SOEMax]))[:,1], Array(normalized_rhs.(model[:SOEMax]))) # adding extra column for T[1]-1
    SOEMin_value = hcat(Array(normalized_rhs.(model[:SOEMin]))[:,1], Array(normalized_rhs.(model[:SOEMin])))
    
    SOEPUP_value = min.(SOEPUP_value, SOEMax_value)
    SOEPDN_value = max.(SOEPDN_value, SOEMin_value)
    return Containers.DenseAxisArray(SOEPUP_value, S, T_incr), Containers.DenseAxisArray(SOEPDN_value, S, T_incr)
    
end

function constrain_SOE_to_envelopes(model, SOEUP_value, SOEDN_value)
    println("Constraining SOE...")
    SOE = model[:SOE]
    S = axes(SOE)[1]
    T_incr = axes(SOE)[2]
    @constraint(model, SOEEnvelopeUP[s in S, t in T_incr],
        SOE[s,t] <= SOEUP_value[s,t]
    )
    @constraint(model, SOEEnvelopeDN[s in S, t in T_incr],
        SOE[s,t] >= SOEDN_value[s,t]
    )
end

function constraint_SOE_final_to_envelopes(model)
    # it assumes envelopes have been calculateed by this stage
    println("Constraining SOE final to envelopes...")
    SOE = model[:SOE]
    T_incr = axes(SOE)[2]
    S = axes(SOE)[1]
    remove_variable_constraint(model, :SOEFinal)
    @constraint(model, SOEFinalDn[s in S],
        SOE[s,T_incr[end]] >= model[:SOEDN_EC][s,T_incr[end]]
    )
    @constraint(model, SOEFinalUp[s in S],
        SOE[s,T_incr[end]] <= model[:SOEUP_EC][s,T_incr[end]]
    )
end

function fix_decision_variables(model, variables, remove_variables_from_objective = true)
    println("Fixing decision variables...")
    for (var, var_value) in variables
    # for (var, var_value) in [(model[var], value.(model[var])) for var in decision_variables]
        for key in collect(keys(var))
           fix(var[key], var_value[key]; force = !is_binary(var[key]))
           if remove_variables_from_objective 
                println("Removing fixed decision variables from objective...")
                set_objective_coefficient(model, var[key], 0) # when set to zero they are removed from objective function
            end 
        end
    end
end

function remove_energy_and_reserve_constraints(model)
    println("Removing reserve, energy reserve and envelope constraints...")
    # Remove reserve, energy reserve and storge envelope's associated variables/constraints
    # TODO: check fix decision variables
    keys = [:ResUpRequirement, :ResDnRequirement,
        :ResUpThermal, :ResDnThermal, :ResUpRamp, :ResDnRamp, :ResUpRampRobust, :ResDnRampRobust,
        :ResUpStorage, :ResDownStorage, 
        :ResUpStorageDisCapacityMax, :ResUpStorageDisLogic, :ResUpStorageChCapacityMax, :ResUpStorageChLogic,
        :ResDownStorageChCapacityMax, :ResDownStorageChLogic, :ResDownStorageDisCapacityMax, :ResDownStorageDisLogic,
        :ResUpStorageDisMax, :ResDownStorageChMax,
        :ResUpStorageCapacityMax, :ResDownStorageCapacityMax,
        # :D,:U,
        :SOEUpEvol, :SOEDnEvol, :SOEUP_0, :SOEDN_0, :SOEUPMax, :SOEDNMax, :SOEUPMin, :SOEDNMin, :SOEDN, :SOEUP,
        # :EnergyResUpThermal, :EnergyResDownThermal, :EnergyResUpRamp, :EnergyResDnRamp,
        # :EnergyResUpZero, :EnergyResDnZero,
        # :EnergyResUpStorage, :EnergyResDownStorage, :EnergyResUpLink, :EnergyResDownLink,
        # :EnergyResUpRequirement, :EnerResDnRequirement,
        :OVMax, :OVMin, :ResUpThermalMin, :ResDownThermalMin, :CommitmentMin, :ResUpStorageMax, :ResUpStorageMin,
        :Startup, :Shutdown, :CommitmentStatus,
        #:RESUP, :RESDN, -> do not remove because belong to the OF.
        :RESUPCH, :RESDNCH, :RESUPDIS, :RESDNDIS,
        # :COMMIT, :START, :SHUT, -> do not remove because belong to the OF.
        ]
    for k in keys
        if haskey(model, k)
            remove_variable_constraint(model, k)
        end
    end
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
    if !is_solved_and_feasible(ed)
        @infiltrate
    end
    if get(kwargs, :save_constraints_status, false) # deprecated
        save_constraints_status(ed, string(get(kwargs, :save_constraints_status_for_demand, nothing)))
    end
    println("done")
    return get_model_solution(ed, gen_df, loads, gen_variable; kwargs...)
end

function solve_economic_dispatch_get_solution(uc, gen_df, loads, gen_variable; kwargs...)
    # Parsing arguments...
    remove_reserve_constraints = get(kwargs, :remove_reserve_constraints, true)
    max_iterations = get(kwargs, :max_iterations, NB_ITERATIONS)
    constrain_dispatch = get(kwargs, :constrain_dispatch, true)
    variables_to_constrain = get(kwargs, :variables_to_constrain, [GEN, CH, DIS])
    remove_variables_from_objective = get(kwargs, :remove_variables_from_objective, false)
    VLOL = get(kwargs, :VLOL, 1e4)
    VLGEN = get(kwargs, :VLGEN, 0)
    reference_solution = get(kwargs, :reference_solution, nothing)
    bidirectional_storage_reserve = get(kwargs, :bidirectional_storage_reserve, true)
    constrain_SOE_by_envelopes = get(kwargs, :constrain_SOE_by_envelopes, false)
    constrain_dispatch_by_energy = get(kwargs, :constrain_dispatch_by_energy, false)
    constraint_dispatch_by_multipliers = get(kwargs, :constraint_dispatch_by_multipliers, false)
    # parsing end
    # uc = construct_unit_commitment(gen_df, loads[!,[HOUR, DEMAND]], gen_variable; kwargs...)
    
    if !isnothing(reference_solution)
        uc = generate_alternative_model(uc, reference_solution)
    end
    # optimize!(uc)
    if constrain_SOE_by_envelopes
        variables_to_constrain = [GEN]
    end
    ed = construct_economic_dispatch(uc, loads[!,[HOUR, DEMAND]], constraint_dispatch_by_multipliers, constrain_SOE_by_envelopes, remove_reserve_constraints, constrain_dispatch, constrain_dispatch_by_energy, bidirectional_storage_reserve, remove_variables_from_objective, variables_to_constrain, VLOL, VLGEN)
    # save_model_to_file(ed,"ed")
    solutions = Dict()
    
    kwargs = Dict(kwargs)
    for k in first(propertynames(loads[!, Not(HOUR)]), max_iterations)
        println("")
        println("Montecarlo iteration: $(k)")
        gen_df_k, loads_df_k, gen_variable_k = pre_process_load_gen_variable(gen_df, rename(loads[!,[HOUR,k]], k=>DEMAND), gen_variable)
        update_demand(ed, loads_df_k)
        update_generation(ed, gen_variable_k)
        if (k == get(kwargs, :save_constraints_status_for_demand, false)) 
            kwargs[:save_constraints_status] = true
        else
            kwargs[:save_constraints_status] = false
        end
        solutions[k] = solve_economic_dispatch_(ed, gen_df_k, loads_df_k, gen_variable_k; kwargs...)
    end
    return merge_solutions(solutions)
end
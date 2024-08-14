using DataFrames

FIELD_FOR_ENRICHING = [:r_id, :resource, :full_id]

function value_to_df(var)
    if var isa JuMP.Containers.DenseAxisArray
        if length(value.(var).axes) == 2
            return value_to_df_2dim(var)
        else
            return value_to_df_1dim(var)
        end
    elseif var isa JuMP.Containers.SparseAxisArray
        return value_to_df_multidim(var, [:r_id, :hour_i, :hour])
    else
        println("Could not identify type of output. Returning initial variable...")
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

function value_to_df_1dim(var)
    return DataFrame(hour = value.(var).axes[1], value = value.(var).data)
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

function get_solution(model)
    return merge(
        get_solution_variables(model),
        (scalar = DataFrame(objective_value = objective_value(model), termination_status = termination_status(model), OPEX = value(model[:OPEX])),)
        # (objective_value = objective_value(model), termination_status = termination_status(model))
    )
end

function get_solution_variables(model)
    variables_to_get = [:GEN, :COMMIT, :SHUT, :CH, :DIS, :SOE, :SOEUP, :SOEDN, :RESUP, :RESUPDIS, :RESUPCH, :RESDNDIS, :RESDNCH, :RESDN, :ERESUP, :ERESDN, :LOL, :LGEN, :SOEUP_EC, :SOEDN_EC]
    return NamedTuple(k => value_to_df(model[k]) for k in intersect(keys(object_dictionary(model)), variables_to_get))
end

function enrich_dfs(solution, gen_df, loads, gen_variable, storage)
    #TODO: deal with missing values
    # Curtailment calculation
    # out = Dict(pairs(solution[[:objective_value, :termination_status]]))
    out = Dict(pairs(solution[[:scalar]]))
    # out = Dict(:generation => get_enriched_generation(solution, gen_df, gen_variable))
    out[:generation] = get_enriched_generation(solution, gen_df, gen_variable)
    out[:generation_parameters] = get_generation_parameters(gen_df)
    
    data = copy(gen_df[!,FIELD_FOR_ENRICHING]) # data for enriching
    if !isnothing(storage)
        append!(data, storage[!,FIELD_FOR_ENRICHING] )
        out[:storage] = get_enriched_storage(solution, data)
        out[:storage_parameters] = get_storage_parameters(storage)
    end
    if haskey(solution, :RESUP) & haskey(solution, :RESDN)
        out[:reserve] =  get_enriched_reserve(solution, data)
    end
    if haskey(solution, :ERESUP) & haskey(solution, :ERESDN)
        out[:energy_reserve] =  get_enriched_energy_reserve(solution, data)
    end
    out[:demand] = get_enriched_demand(solution, loads)
    return NamedTuple(out)
end


function get_enriched_energy_reserve(solution, data)
    return leftjoin(
        innerjoin(
            rename(solution.ERESUP, :value => :reserve_up_MW),
            rename(solution.ERESDN, :value => :reserve_down_MW),
            on = [:r_id, :hour, :hour_i]
        ),
        data[!,FIELD_FOR_ENRICHING],
        on = :r_id
    )
end

function get_enriched_reserve(solution, data)
    aux = leftjoin(
        innerjoin(
            rename(solution.RESUP, :value => :reserve_up_MW),
            rename(solution.RESDN, :value => :reserve_down_MW),
            on = [:r_id, :hour]
        ),
        data[!,FIELD_FOR_ENRICHING],
        on = :r_id
    )
    if haskey(solution, :RESUPDIS)
        aux = outerjoin(
            aux,
            rename(solution.RESUPDIS, :value => :reserve_discharge_up_MW),
            rename(solution.RESUPCH, :value => :reserve_charge_up_MW),
            rename(solution.RESDNDIS, :value => :reserve_discharge_down_MW),
            rename(solution.RESDNCH, :value => :reserve_charge_down_MW),
            on = [:r_id, :hour],
        )
    end
    return aux
end

function get_enriched_storage(solution, data)
    aux = innerjoin(
        rename(solution.CH, :value => :charge_MW),
        rename(solution.DIS, :value => :discharge_MW),
        rename(solution.SOE, :value => :SOE_MWh),
        on = [:r_id, :hour]
    )
    if haskey(solution, :SOEUP) & haskey(solution, :SOEDN)
        aux = innerjoin(
            aux,
            rename(solution.SOEUP, :value => :envelope_up_MWh),
            rename(solution.SOEDN, :value => :envelope_down_MWh),
            on = [:r_id, :hour]
        )
    end
    if haskey(solution, :SOEUP_EC) & haskey(solution, :SOEDN_EC)
        aux = innerjoin(
            aux,
            rename(solution.SOEUP_EC, :value => :envelope_up_MWh),
            rename(solution.SOEDN_EC, :value => :envelope_down_MWh),
            on = [:r_id, :hour]
        )
    end
    return leftjoin(aux, data[!,FIELD_FOR_ENRICHING], on = :r_id)
end

function get_enriched_generation(solution, gen_df, gen_variable)
    curtail = innerjoin(gen_variable, solution.GEN, on = [:r_id, :hour])
    curtail.value = curtail.cf .* curtail.existing_cap_mw - curtail.value
    aux = outerjoin(
        outerjoin(
            rename(solution.GEN, :value => :production_MW),
            rename(curtail[!,[:r_id, :hour, :value]], :value => :curtailment_MW),
            rename(solution.COMMIT, :value => :commit),
            on = [:r_id, :hour]
        ),
        gen_df[!,FIELD_FOR_ENRICHING],
        on = :r_id
    )
    replace!(aux.curtailment_MW, missing => 0)
    return aux
end

function get_enriched_demand(solution, loads)
    demand = rename(loads, [ x => "demand_MW" for x in names(loads[!,Not(:hour)])])
    # demand =  rename(loads, :demand => :demand_MW)
    demand.r_id .= missing
    demand.resource .= "total"
    if haskey(solution, :LOL)
        demand.demand_MW =  demand.demand_MW - solution[:LOL].value
        demand = leftjoin(
            demand, 
            rename(solution[:LOL], :value => :LOL_MW),
            on = [:hour])
    end
    if haskey(solution, :LGEN)
        LGEN = copy(solution.LGEN)
        # LGEN.r_id .= missing
        LGEN.resource .= "total"
        demand = outerjoin(
            demand,
            rename(LGEN[!,[:resource, :hour, :value]], :value => :LGEN_MW),
            on = [:hour, :resource]
        )
        replace!(demand.LGEN_MW, missing => 0)
    end
    return demand
end

function get_generation_parameters(gen_df)
    parameters_to_get = [:existing_cap_mw, :min_power]
    return rename(copy(gen_df[!,union(FIELD_FOR_ENRICHING, parameters_to_get)]), :existing_cap_mw => :P_max_MW)
end

function get_storage_parameters(storage)
    parameters_to_get = [:existing_cap_mw, :max_energy_mwh, :charge_efficiency, :discharge_efficiency]
    return rename(storage[!,union(FIELD_FOR_ENRICHING, parameters_to_get)],[:existing_cap_mw, :max_energy_mwh] .=> [:P_max_MW, :SOE_max_MWh])
end

function get_model_solution(model, gen_df, loads, gen_variable; config...)
    # Model is either UC or ED
    if get(config, :enriched_solution, true)
        return enrich_dfs(get_solution(model), gen_df, loads, gen_variable, get(config, :storage, nothing)) 
    else
        return get_solution(model)#TODO: remove kwargs 
    end
end
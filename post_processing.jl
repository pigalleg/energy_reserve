using DataFrames

FIELD_FOR_ENRICHING = [:r_id, :resource, :full_id]

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

function get_solution_variables(model)
    variables_to_get = [:GEN, :COMMIT, :SHUT, :CH, :DIS, :SOE, :RESUP, :RESDN, :ERESUP, :ERESDN]
    return NamedTuple(k => value_to_df(model[k]) for k in intersect(keys(object_dictionary(model)), variables_to_get))
end

# function add_reserve_to_solution(enriched_solution_value, solution)
#     reserve = innerjoin(
#         rename(solution.RESUP, :value => :reserve_up_MW),
#         rename(solution.RESDN, :value => :reserve_down_MW),
#         on = [:r_id, :hour]
#     )
#     # leftjoin allows to add reserve to frame in the left
#     return leftjoin(enriched_solution_value, reserve, on = [:r_id, :hour])
# end

# function add_energy_reserve_to_solution(enriched_solution_value, solution)
#     reserve = innerjoin(
#         rename(solution.ERESUP, :value => :reserve_up_MW),
#         rename(solution.ERESDN, :value => :reserve_down_MW),
#         on = [:r_id, :hour, :hour_i]
#     )
#     return outerjoin(enriched_solution_value, reserve, on = [:r_id, :hour])
# end
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
    return leftjoin(
        innerjoin(
            rename(solution.RESUP, :value => :reserve_up_MW),
            rename(solution.RESDN, :value => :reserve_down_MW),
            on = [:r_id, :hour]
        ),
        data[!,FIELD_FOR_ENRICHING],
        on = :r_id
    )
end

function get_enriched_storage(solution, data)
    aux = 
    return leftjoin(
        innerjoin(
            rename(solution.CH, :value => :charge_MW),
            rename(solution.DIS, :value => :discharge_MW),
            rename(solution.SOE, :value => :SOE_MWh),
            on = [:r_id, :hour]
        ),
        data[!,FIELD_FOR_ENRICHING],
        on = :r_id
    )
end

function get_enriched_generation(solution, gen_df, gen_variable)
    curtail = innerjoin(gen_variable, solution.GEN, on = [:r_id, :hour])
    curtail.value = curtail.cf .* curtail.existing_cap_mw - curtail.value
    aux = outerjoin(
        outerjoin(
            rename(solution.GEN, :value => :production_MW),
            rename(curtail[!,[:r_id, :hour, :value]], :value => :curtailment_MW), 
            on = [:r_id, :hour]
        ),
        gen_df[!,FIELD_FOR_ENRICHING],
        on = :r_id
    )
    replace!(aux.curtailment_MW, missing => 0)
    return aux
end

function get_enriched_demand(loads)
    demand =  rename(loads, :demand => :demand_MW)
    demand.r_id .= missing
    demand.resource .= "total"
    return demand
end

function to_enriched_df(solution, gen_df, loads, gen_variable; kwargs...)
    #TODO: deal with missing values
    # Curtailment calculation
    out = Dict(:generation => get_enriched_generation(solution, gen_df, gen_variable))
    data = copy(gen_df[!,FIELD_FOR_ENRICHING]) # data for enriching
    if haskey(kwargs, :storage)
        append!(data, kwargs[:storage][!,FIELD_FOR_ENRICHING] )
        out[:storage] = get_enriched_storage(solution, data)
    end
    if haskey(solution, :RESUP) & haskey(solution, :RESDN)
        out[:reserve] =  get_enriched_reserve(solution, data)
    end
    if haskey(solution, :ERESUP) & haskey(solution, :ERESDN)
        out[:energy_reserve] =  get_enriched_energy_reserve(solution, data)
    end
    out[:demand] = get_enriched_demand(loads)
    return NamedTuple(out)
end 
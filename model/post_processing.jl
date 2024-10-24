using DataFrames

FIELD_FOR_ENRICHING = [:r_id, :resource, :full_id]
SOLUTION_KEYS = [:demand, :generation, :storage, :reserve, :energy_reserve, :scalar, :generation_parameters, :storage_parameters ,:objective_function]

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
    variables_to_get = [:GEN, :COMMIT, :SHUT, :START, :CH, :DIS, :SOE, :SOEUP, :SOEDN, :RESUP, :CHRESDNCH, :CHRESUPCH, :DISRESUPDIS, :DISRESDNDIS, :RESDN, :ERESUP, :ERESDN, :LOL, :LGEN, :SOEUP_EC, :SOEDN_EC]
    return NamedTuple(k => value_to_df(model[k]) for k in intersect(keys(object_dictionary(model)), variables_to_get))
end

function enrich_dfs(solution, gen_df, loads, gen_variable, storage, parameters)
    #TODO: deal with missing values
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
    out[:objective_function] = get_enriched_objective_value(out, gen_df, storage,parameters)
    # out[:constraints] = get_constraints(solution)
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
    # if haskey(solution, :RESUPDIS)
    #     aux = outerjoin(
    #         aux,
    #         rename(solution.RESUPDIS, :value => :reserve_discharge_up_MW),
    #         rename(solution.RESUPCH, :value => :reserve_charge_up_MW),
    #         rename(solution.RESDNDIS, :value => :reserve_discharge_down_MW),
    #         rename(solution.RESDNCH, :value => :reserve_charge_down_MW),
    #         on = [:r_id, :hour],
    #     )
    # end
    if haskey(solution, :CHRESDNCH)
        aux = outerjoin(
            aux,
            rename(solution.CHRESDNCH, :value => :chargue_reserve_max_MW),
            rename(solution.CHRESUPCH, :value => :chargue_reserve_min_MW),
            rename(solution.DISRESUPDIS, :value => :discharge_reserve_max_MW),
            rename(solution.DISRESDNDIS, :value => :discharge_reserve_min_MW),
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
    # if haskey(solution, :SOEUP_ED) & haskey(solution, :SOEDN_EC) # deprecated
    #     aux = innerjoin(
    #         aux,
    #         rename(solution.SOEUP_ED, :value => :envelope_up_MWh),
    #         rename(solution.SOEDN_ED, :value => :envelope_down_MWh),
    #         on = [:r_id, :hour]
    #     )
    # end
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
            rename(solution.START, :value => :start),
            rename(solution.SHUT, :value => :shut),
            on = [:r_id, :hour]
        ),
        gen_df[!,FIELD_FOR_ENRICHING],
        on = :r_id
    )
    replace!(aux.curtailment_MW, missing => 0)
    return aux
end

function get_enriched_demand(solution, loads)
    demand_ = "day" in names(loads) ? select(loads,Not(:day)) : loads
    demand = rename(demand_, [ x => "demand_MW" for x in names(select(demand_,Not(:hour)))])
    # demand =  rename(loads, :demand => :demand_MW)
    demand.r_id .= missing
    demand.resource .= "system"
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
        LGEN.resource .= "system"
        demand = outerjoin(
            demand,
            rename(LGEN[!,[:resource, :hour, :value]], :value => :LGEN_MW),
            on = [:hour, :resource]
        )
        replace!(demand.LGEN_MW, missing => 0)
    end
    return demand
end

function get_enriched_objective_value(enriched_solution, gen_df, storage, parameters)
    function check()
        within_MIPGap(x,y,tolerance) = abs((x-y)/x) <= tolerance
        sum_cost = sum(skipmissing(cost.production_cost)) + sum(skipmissing(cost.start_cost))
        if !within_MIPGap(enriched_solution[:scalar].OPEX[1], sum_cost , parameters.MIPGap) 
            error("Start and operational cost missmatch with OPEX")
        end
        if :reserve_cost in propertynames(cost)
            sum_cost += sum(skipmissing(cost.reserve_cost))
        end
        if :LOL_cost in propertynames(cost)
            sum_cost += sum(skipmissing(cost.LOL_cost))  + sum(skipmissing(cost.LGEN_cost))  
        end
        if !within_MIPGap(enriched_solution[:scalar].objective_value[1], sum_cost, parameters.MIPGap) 
            error("Start, operational cost and reserve penalization missmatch with objective value")
        end
    end
   
    cost_fields = [:start_cost_per_mw, :existing_cap_mw, :heat_rate_mmbtu_per_mwh, :fuel_cost, :var_om_cost_per_mwh, :fixed_om_cost_per_mw_per_hour]
    fields_to_remove = [:production_MW, :curtailment_MW, :commit, :start, :full_id, :shut]
    cost = leftjoin(
        enriched_solution[:generation],
        gen_df[!,union(cost_fields, [:r_id])],
        on = [:r_id]
        )

    cost.start_cost = cost.start_cost_per_mw .* cost.existing_cap_mw .* cost.start
    cost.production_cost = 
        (cost.heat_rate_mmbtu_per_mwh .* cost.fuel_cost + cost.var_om_cost_per_mwh) .* cost.production_MW +
        cost.fixed_om_cost_per_mw_per_hour .* cost.existing_cap_mw .* replace(cost.commit, missing => 1)
    select!(cost,Not(union(cost_fields,fields_to_remove)))
    
    if :storage in keys(enriched_solution)
        cost_fields = [:var_om_cost_per_mwh]
        fields_to_remove = intersect([:charge_MW, :discharge_MW, :SOE_MWh, :envelope_up_MWh,:envelope_down_MWh, :full_id], propertynames(enriched_solution[:storage]))  
        storage_cost = leftjoin(
            enriched_solution[:storage],
            storage[!,union(cost_fields, [:r_id])],
            on = [:r_id]
        )
        storage_cost.production_cost =  (storage_cost.charge_MW + storage_cost.discharge_MW) .* storage_cost.var_om_cost_per_mwh
        select!(storage_cost, Not(union(cost_fields,fields_to_remove)))
        cost = vcat(cost, storage_cost, cols=:union)
    end
    if :reserve in keys(enriched_solution)
        fields_to_remove = [:reserve_up_MW, :reserve_down_MW, :full_id]
        reserve_cost = copy(enriched_solution[:reserve])
        reserve_cost.reserve_cost = (reserve_cost.reserve_up_MW + reserve_cost.reserve_down_MW)*parameters.VRESERVE
        select!(reserve_cost, Not(fields_to_remove)) 
        cost = vcat(cost, reserve_cost, cols=:union)
    end

    if :energy_reserve in keys(enriched_solution)
        # This one gets reserve cost for energy reserve on the diagonal terms!
        fields_to_remove = [:reserve_up_MW, :reserve_down_MW, :full_id, :hour_i]
        reserve_cost = filter(y->(y.hour_i .== y.hour),  enriched_solution[:energy_reserve])
        reserve_cost.reserve_cost = (reserve_cost.reserve_up_MW + reserve_cost.reserve_down_MW)*parameters.VRESERVE
        select!(reserve_cost, Not(fields_to_remove)) 
        cost = vcat(cost, reserve_cost, cols=:union)
    end
    if :LOL_MW in propertynames(enriched_solution[:demand]) # ED, we assume that LGEN_MW is also present when LOL_MW is present
        fields_to_remove = [:LOL_MW, :LGEN_MW, :demand_MW]
        losses_cost = copy(enriched_solution[:demand])
        losses_cost.LOL_cost = losses_cost.LOL_MW .* parameters.VLOL
        losses_cost.LGEN_cost = losses_cost.LGEN_MW .* parameters.VLGEN
        select!(losses_cost, Not(fields_to_remove))
        cost = vcat(cost, losses_cost, cols=:union) 
    end
    check()
    return cost
end

function get_generation_parameters(gen_df)
    parameters_to_get = [:existing_cap_mw, :min_power]
    return rename(copy(gen_df[!,union(FIELD_FOR_ENRICHING, parameters_to_get)]), :existing_cap_mw => :P_max_MW)
end

function get_storage_parameters(storage)
    parameters_to_get = [:existing_cap_mw, :max_energy_mwh, :charge_efficiency, :discharge_efficiency, :initial_energy_proportion]
    return rename(storage[!,union(FIELD_FOR_ENRICHING, parameters_to_get)],[:existing_cap_mw, :max_energy_mwh] .=> [:P_max_MW, :SOE_max_MWh])
end

function get_model_solution(model, gen_df, loads, gen_variable; config...)
    # Model is either UC or ED
    parameters_for_enriching = (MIPGap = parameter_value(model[:MIPGap]), VRESERVE = parameter_value(model[:VRESERVE]))
    if haskey(model, :VLOL) # ED
        parameters_for_enriching = merge(parameters_for_enriching,
            (VLOL = Array(parameter_value.(model[:VLOL])),
            VLGEN = Array(parameter_value.(model[:VLGEN])))
        )
    end
    if get(config, :enriched_solution, true)
        return enrich_dfs(get_solution(model), gen_df, loads, gen_variable, get(config, :storage, nothing),parameters_for_enriching) 
    else
        return get_solution(model)#TODO: remove kwargs 
    end
end

function solution_to_parquet(s, file_name, file_folder)
    # TODO move to post_processing
    if !isdir(file_folder) mkdir(file_folder) end
    println("writing...")
    for (k,v) in zip(propertynames(s), s)
      println("$(file_name)_$k")
      Parquet2.writefile(joinpath(file_folder, file_name*"_"*string(k)*".parquet"), change_type(change_type(v, Symbol, string), TerminationStatusCode, string))
    end
    println("...done")
  end
  
  function parquet_to_solution(file_name, file_folder)
    # TODO 1 convert to TerminationStatusCode
    # TODO 2 move to post_processing
    keys = [k for k in SOLUTION_KEYS if isfile(joinpath(file_folder, file_name*"_"*string(k)*".parquet"))]
    println("reading...")
    aux = [read_parquet_and_convert(joinpath(file_folder, file_name*"_"*string(k)*".parquet")) for k in keys]
    println("...done")
    return NamedTuple(keys .=> aux)
  end
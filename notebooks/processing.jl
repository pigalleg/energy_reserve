using DataFrames

function generate_reserves(loads, margin_percentage, baseload = 0)
  required_reserve = DataFrame(
    hour = loads[!,:hour],
    reserve_up_MW = baseload .+ loads[!,:demand].*margin_percentage,
    reserve_down_MW = loads[!,:demand].*margin_percentage)

  required_energy_reserve = [(row_1.hour, row_2.hour, row_1.reserve_up_MW*(row_1.hour == row_2.hour), row_1.reserve_down_MW*(row_1.hour == row_2.hour)) for row_1 in eachrow(required_reserve), row_2 in eachrow(required_reserve) if row_1.hour <= row_2.hour]
  required_energy_reserve = DataFrame(required_energy_reserve)
  required_energy_reserve = rename(required_energy_reserve, :1 => :i_hour, :2 => :t_hour, :3 => :reserve_up_MW, :4 => :reserve_down_MW,)

  required_energy_reserve_cumulated = [(row_1.hour, row_2.hour, sum(required_reserve[(required_reserve.hour .>= row_1.hour).&(required_reserve.hour .<= row_2.hour),:reserve_up_MW]), sum(required_reserve[(required_reserve.hour .>= row_1.hour).&(required_reserve.hour .<= row_2.hour),:reserve_down_MW])) for row_1 in eachrow(required_reserve), row_2 in eachrow(required_reserve) if row_1.hour <= row_2.hour]
  required_energy_reserve_cumulated = DataFrame(required_energy_reserve_cumulated)
  required_energy_reserve_cumulated = rename(required_energy_reserve_cumulated, :1 => :i_hour, :2 => :t_hour, :3 => :reserve_up_MW, :4 => :reserve_down_MW,)

  return required_reserve, required_energy_reserve, required_energy_reserve_cumulated
end


function calculate_supply_demand(solution, group_by = [:hour, :resource] )
  #Supply-demand computation
  # group_by = intersect(propertynames(solution.generation),[:hour, :resource, :iteration])

  supply = combine(groupby(solution.generation, group_by), :production_MW => sum, renamecols=false)
  demand = combine(groupby(solution.demand, group_by), :demand_MW => sum, renamecols=false)


  aux = combine(groupby(solution.generation, group_by), :curtailment_MW => sum, renamecols=false)
  # replace!(aux.curtailment_MW, missing => 0)
  aux = aux[aux.curtailment_MW.>0,:]
  rename!(aux, :curtailment_MW => :production_MW)
  transform!(aux, :resource .=> ByRow(x -> x*"_curtailment") => :resource)
  append!(supply, aux, promote = true)

  if haskey(solution,:storage)
      aux = combine(groupby(solution.storage, group_by), [:discharge_MW => sum, :charge_MW => sum], renamecols=false)
      rename!(aux, [:discharge_MW => :production_MW, :charge_MW => :demand_MW])
      append!(supply, aux[!, push!(copy(group_by), :production_MW)])
      append!(demand,  aux[!,push!(copy(group_by), :demand_MW)], promote = true)
  end 
  return supply, demand
end

function calculate_reserve(reserve, required_reserve, group_by_ = [:hour, :resource])
  field_up_dn = [:reserve_up_MW,:reserve_down_MW]
  aux = reserve[!,union(group_by_, field_up_dn)]
  replace!([aux.reserve_up_MW, missing => 0, aux.reserve_up_MW, missing => 0])
  # replace!(aux.reserve_down_MW, missing => 0)
  group_by = intersect(propertynames(aux), group_by_)
  reserve = combine(groupby(aux, group_by), [field_up_dn[1] => sum, field_up_dn[2] => sum], renamecols=false)
   # TODO: Adapt the following lines to consider the cases where group_by has more keys (e.g. :iteration, :coniguration)
  aux = required_reserve
  aux.resource.= "required"
  append!(reserve, aux)
  return reserve
end

function calculate_battery_reserve(solution_storage, solution_reserve, efficiency = 0.92)
  # TODO: adapt when input dfs have more keys (e.g., :iteration, configuration)
  variables_to_get = [:r_id, :hour, :SOE_MWh, :reserve_up_MW, :reserve_down_MW, :envelope_up_MWh, :envelope_down_MWh, :iteration]
  
  fields_storage = intersect(propertynames(solution_storage), variables_to_get)
  fields_solution_reserve = intersect(propertynames(solution_reserve), variables_to_get)
  group_by = intersect([:r_id, :hour, :iteration], fields_storage, fields_solution_reserve)
  out = innerjoin(
      solution_storage[!,fields_storage], 
      solution_reserve[solution_reserve.resource .=="battery", fields_solution_reserve], 
      on = intersect(fields_storage, fields_storage, group_by))
  out = combine(groupby(out, setdiff(group_by, [:r_id])), Not(setdiff(group_by, [:r_id])) .=> sum, renamecols=false)
  out.reserve_down_MW_eff .= out.reserve_down_MW*efficiency
  out.reserve_up_MW_eff .= -out.reserve_up_MW*1/efficiency
  return out
end

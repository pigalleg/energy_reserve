using DataFrames
using Parquet2
using MathOptInterface: TerminationStatusCode
order = [
  "solar_photovoltaic_curtailment",
  "onshore_wind_turbine_curtailment",
  "small_hydroelectric_curtailment",
  "loss_of_generation",
  "total_loss_of_load_ED",
  "total_loss_of_generation_ED",
  "net_generation_curtailment",
  "battery",
  "solar_photovoltaic",
  "net_generation",
  "natural_gas_fired_combustion_turbine",
  "natural_gas_fired_combined_cycle",
  "onshore_wind_turbine",
  "hydroelectric_pumped_storage",
  "small_hydroelectric",
  "biomass",
  "total",
  "required"]


G_NET_GENERAION_FULL_ID = "net_generation"

function order_df(df_)
  df = copy(df_)
  df[!, :order] = indexin(df[!,:resource], order)
  replace!(df[!,:order], nothing =>length(order))
  sort!(df, :order, rev = true)
  return select(df, Not(:order))
end

function generate_reserves(loads, gen_variable, margin_percentage, baseload = 0)
  filter = gen_variable[!,:full_id] .== G_NET_GENERAION_FULL_ID
  net_gen = gen_variable[filter,:cf] .* gen_variable[filter,:existing_cap_mw]
  
  required_reserve = DataFrame(
    hour = loads[!,:hour],
    reserve_up_MW = baseload .+ (loads[!,:demand] .+ net_gen).*margin_percentage,
    reserve_down_MW = (loads[!,:demand] .+ net_gen).*margin_percentage)
  required_reserve = (required_reserve.>0).*required_reserve .- (required_reserve.<0).*required_reserve # negative to positive values

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
  demand = combine(groupby(solution.demand, group_by), :demand_MW => sum, renamecols=false)
  # replace!(aux.curtailment_MW, missing => 0)
  if :LOL_MW in propertynames(solution.demand)
    aux = combine(groupby(solution.demand, group_by), :LOL_MW => sum, renamecols=false)
    aux = aux[aux.LOL_MW.>0,:]
    rename!(aux, :LOL_MW => :demand_MW)
    transform!(aux, :resource .=> ByRow(x -> x*"_loss_of_load_ED") => :resource)
    append!(demand, aux, promote = true)
  end
  if :LGEN_MW in propertynames(solution.demand)
    aux = combine(groupby(solution.demand, group_by), :LGEN_MW => sum, renamecols=false)
    aux = aux[aux.LGEN_MW.>0,:]
    rename!(aux, :LGEN_MW => :demand_MW)
    transform!(aux, :resource .=> ByRow(x -> x*"_loss_of_generation_ED") => :resource)
    append!(demand, aux, promote = true)
  end
  supply = combine(groupby(solution.generation, group_by), :production_MW => sum, renamecols=false)
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
  return order_df(supply), order_df(demand)
end

function calculate_reserve(reserve, required_reserve = nothing, group_by_ = [:hour, :resource])
  field_up_dn = [:reserve_up_MW,:reserve_down_MW]
  aux = reserve[!,union(group_by_, field_up_dn)]
  replace!([aux.reserve_up_MW, missing => 0, aux.reserve_up_MW, missing => 0])
  # replace!(aux.reserve_down_MW, missing => 0)
  group_by = intersect(propertynames(aux), group_by_)
  reserve = combine(groupby(aux, group_by), [field_up_dn[1] => sum, field_up_dn[2] => sum], renamecols=false)
   # TODO: Adapt the following lines to consider the cases where group_by has more keys (e.g. :iteration, :coniguration)
  if !isnothing(required_reserve) 
    aux = copy(required_reserve)
    aux.resource.= "required"
    append!(reserve, aux)
  end
  return order_df(reserve)
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

function change_type(df, from, to)
  return mapcols(x -> eltype(x) == from ? to.(x) : x, df)
end


parse_configuration_to_mu(x) = !isnothing(match(r"base_ramp_storage_envelopes_up_(\w+)_dn_(\w+)", string(x))) ? parse(Float64, replace(match(r"base_ramp_storage_envelopes_up_(\w+)_dn_(\w+)", string(x))[1], "_" => ".")) : missing

function calculate_adecuacy_gcdi_KPI(s_ed, s_uc, thres =.001) # thres = 1 Watt
  f_LOL(x,y) = 
    (LLD_h=count(x.>thres),
    # LOLP=count(x.>thres)/length(y)*100,
    ENS_MWh = sum(x),
    # LOL_percentage = sum(x)/(sum(y) + sum(x))*100,
    demand_MWh = (sum(y) + sum(x)),
    )

  f_CUR(x,y) =     
      (CURD_h=count(x.>thres),
      # CURP=count(x.>thres)/length(y)*100,
      CUR_MWh = sum(x),
      # CUR_percentage = sum(x)/(sum(y) + sum(x))*100,
      RES_production_MWh = (sum(y) + sum(x)),
      )

  group_by_big = [:configuration, :day, :iteration]
  filter = in(["onshore_wind_turbine","small_hydroelectric","solar_photovoltaic", "net_generation"]).(s_ed.generation.resource)
  
  gcdi_KPI = outerjoin(
    combine(groupby(s_ed.demand, group_by_big), [:LOL_MW, :demand_MW] => ((x,y)->f_LOL(x,y)) => AsTable), 
    combine(groupby(s_ed.generation[filter,:], group_by_big), [:curtailment_MW, :production_MW] =>((x,y) -> f_CUR(x,y))=> AsTable),
    # combine(groupby(s_ed.demand, group_by_big), :LGEN_MW => sum => :LGEN_MWh),
    on = group_by_big)
  if :LGEN_MW in propertynames(s_ed.demand)
    gcdi_KPI = outerjoin(gcdi_KPI, combine(groupby(s_ed.demand, group_by_big), :LGEN_MW => sum => :LGEN_MWh),on = group_by_big)
  end
  # folowing leftjoin will repeat values right values for "iteration"
  s_ed_scalar = s_ed.scalar[:, Not(:termination_status)]
  leftjoin!(
      s_ed_scalar, 
      rename(s_uc.scalar[:, Not(:termination_status)], :objective_value =>:objective_value_uc),
      on = setdiff(propertynames(s_ed_scalar), [:objective_value, :iteration]))
      
  # relative difference with respect to s_uc
  s_ed_scalar.Δobjective_value = s_ed_scalar.objective_value .- s_ed_scalar.objective_value_uc
  s_ed_scalar.Δobjective_value_relative = (s_ed_scalar.objective_value .- s_ed_scalar.objective_value_uc)./s_ed_scalar.objective_value_uc

  # relative difference with respect to ref_configuration = :base_ramp_storage_envelopes_up_0_dn_0
  s_ed_scalar = include_Δobjective_value(s_ed_scalar)

  leftjoin!(gcdi_KPI, s_ed_scalar, on = group_by_big)
  transform!(gcdi_KPI, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu)
  sort!(gcdi_KPI, :mu)
  return gcdi_KPI
end

function include_Δobjective_value(s_scalar, keys = [:day, :iteration], reference_configuration_key = :base_ramp_storage_envelopes_up_0_dn_0)
  ref_configuration = s_scalar[s_scalar.configuration .== reference_configuration_key, union([:objective_value], keys)]

  s_scalar = leftjoin(s_scalar, 
      rename(ref_configuration, :objective_value => :objective_value_ref_conf), 
      on = keys)

  s_scalar.Δobjective_value_ref_conf = s_scalar.objective_value .- s_scalar.objective_value_ref_conf
  s_scalar.Δobjective_value_relative_ref_conf = (s_scalar.objective_value .- s_scalar.objective_value_ref_conf)./s_scalar.objective_value_ref_conf
  return s_scalar
end

function calculate_adecuacy_gcd_KPI(gcdi_KPI)
  group_by = [:configuration, :day]
  gcd_KPI = outerjoin(
    combine(groupby(gcdi_KPI, group_by), [:LLD_h, :ENS_MWh] => ((x,y)->(LOLE = mean(x), EENS = mean(y))) => AsTable),  #TODO: change format
    combine(groupby(gcdi_KPI, group_by), [:CURD_h, :CUR_MWh] => ((x,y)->(CURE = mean(x), ECUR = mean(y))) => AsTable), #TODO: change format
    combine(groupby(gcdi_KPI, group_by), [:objective_value, :Δobjective_value_relative_ref_conf] .=> mean .=> [:EOV, :EΔOV]),
    on=[:configuration, :day])
  if :LGEN_MWh in propertynames(gcdi_KPI)
    gcd_KPI = outerjoin(gcd_KPI, combine(groupby(gcdi_KPI, group_by), :LGEN_MWh => sum => :ELGEN),on = [:configuration, :day])
  end
  transform!(gcd_KPI, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu)
  sort!(gcd_KPI, :mu)
  return gcd_KPI
end

function calculate_reserve_KPI(s_ed, s_uc)
  group_by = [:configuration, :day]
  group_by_big = [:configuration, :day, :iteration]
  KPI_reserve = leftjoin(
    s_ed.demand[!,union(group_by_big,[:hour, :demand_MW, :LOL_MW])], 
    rename(s_uc.demand[!, union(group_by,[:hour, :demand_MW])], :demand_MW => :demand_uc_MW),
    on = union(group_by,[:hour]))

  required = (KPI_reserve.demand_MW .+ KPI_reserve.LOL_MW.- KPI_reserve.demand_uc_MW)
  KPI_reserve.required_r_up_MW = required.*(required.>=0)
  KPI_reserve.required_r_dn_MW = - required.*(required.<0)

  KPI_reserve = outerjoin(
    KPI_reserve,
    combine(groupby(s_ed.generation, union(group_by_big, [:hour])), :production_MW => sum => :production_MW),
    combine(groupby(s_ed.storage, union(group_by_big, [:hour])), [:charge_MW, :discharge_MW] .=> sum .=> [:charge_MW, :discharge_MW]),
    on = union(group_by,[:iteration, :hour]))

  delivered = KPI_reserve.production_MW .+ KPI_reserve.discharge_MW .- KPI_reserve.charge_MW .- KPI_reserve.demand_uc_MW
  KPI_reserve.delivered_r_up_MW = delivered.*(delivered.>=0)
  KPI_reserve.delivered_r_dn_MW = -(delivered .+ KPI_reserve.LOL_MW).*(delivered.<0).*(required.<0)

  KPI_reserve.delivered_r_up_ratio  = KPI_reserve.delivered_r_up_MW./KPI_reserve.required_r_up_MW
  KPI_reserve.delivered_r_dn_ratio  = KPI_reserve.delivered_r_dn_MW./KPI_reserve.required_r_dn_MW
  return KPI_reserve
end

function calculate_reserve_gcdi_KPI(KPI_reserve)
  group_by_big = [:configuration, :day, :iteration]
  gcdi_KPI_reserve  = combine(groupby(KPI_reserve, group_by_big), [:required_r_up_MW, :required_r_dn_MW, :delivered_r_up_MW, :delivered_r_dn_MW] .=> sum .=> [:required_r_up_MWh, :required_r_dn_MWh, :delivered_r_up_MWh, :delivered_r_dn_MWh])
  gcdi_KPI_reserve.delivered_r_up_ratio  = gcdi_KPI_reserve.delivered_r_up_MWh ./gcdi_KPI_reserve.required_r_up_MWh
  gcdi_KPI_reserve.delivered_r_dn_ratio  = gcdi_KPI_reserve.delivered_r_dn_MWh ./gcdi_KPI_reserve.required_r_dn_MWh
  transform!(gcdi_KPI_reserve, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu)
  sort!(gcdi_KPI_reserve, :mu)
  return gcdi_KPI_reserve
end

function calculate_reserve_gcd_KPI(gcdi_KPI_reserve)
  group_by = [:configuration, :day]
  gcd_KPI_reserve  = combine(groupby(gcdi_KPI_reserve, group_by), 
    [:required_r_up_MWh, :required_r_dn_MWh, :delivered_r_up_MWh, :delivered_r_dn_MWh] .=> mean .=> [:E_required_r_up_MWh, :E_required_r_dn_MWh, :E_delivered_r_up_MWh, :E_deliverered_dn_MWh])
  gcd_KPI_reserve.delivered_r_up_ratio  = gcd_KPI_reserve.E_delivered_r_up_MWh ./gcd_KPI_reserve.E_required_r_up_MWh
  gcd_KPI_reserve.delivered_r_dn_ratio  = gcd_KPI_reserve.E_deliverered_dn_MWh ./gcd_KPI_reserve.E_required_r_dn_MWh
  transform!(gcd_KPI_reserve, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu)
  sort!(gcd_KPI_reserve, :mu)
  return gcd_KPI_reserve
end
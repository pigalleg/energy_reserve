using DataFrames
using Parquet2
using MathOptInterface: TerminationStatusCode
using Statistics
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
  "system",
  "required"]


G_NET_GENERAION_FULL_ID = "net_generation"

function order_df(df_)
  df = copy(df_)
  df[!, :order] = indexin(df[!,:resource], order)
  replace!(df[!,:order], nothing =>length(order))
  sort!(df, :order, rev = true)
  return select(df, Not(:order))
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
      aux = combine(groupby(coalesce.(solution.storage,0), group_by), [:discharge_MW => sum, :charge_MW => sum], renamecols=false) #we can do coalesce because we are summing
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


parse_configuration_to_mu(x) = !isnothing(match(r"base_ramp_storage_envelopes_up_(\w+)_dn_(\w+)", string(x))) ? parse(Float64, replace(match(r"base_ramp_storage_envelopes_up_(\w+)_dn_(\w+)", string(x))[1], "_" => ".")) : 1

function calculate_adecuacy_gcdi_KPI(s_ed, s_uc, thres =.001) # thres = 1 Watt
  f_LOL(x,y) = 
    (LLD_h=count(x.>thres),
    # LOLP=count(x.>thres)/length(y)*100,
    ENS_MWh = sum(x),
    # LOL_percentage = sum(x)/(sum(y) + sum(x))*100,
    # demand_MWh = (sum(y) + sum(x)),
    )
  f_CUR(x,y) =     
      (CURD_h=count(x.>thres),
      # CURP=count(x.>thres)/length(y)*100,
      CUR_MWh = sum(x),
      # CUR_percentage = sum(x)/(sum(y) + sum(x))*100,
      # RES_production_MWh = (sum(y) + sum(x)),
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
  # beging REMOVE -- 
  # relative difference with respect to s_uc
  # s_ed_scalar.Δobjective_value = s_ed_scalar.objective_value .- s_ed_scalar.objective_value_uc
  # s_ed_scalar.Δobjective_value_relative = (s_ed_scalar.objective_value .- s_ed_scalar.objective_value_uc)./s_ed_scalar.objective_value_uc

  # relative difference with respect to ref_configuration = :base_ramp_storage_envelopes_up_0_dn_0
  # s_ed_scalar = include_Δobjective_value(s_ed_scalar)
  # end REMOVE -- 

  leftjoin!(gcdi_KPI, calculate_objective_function_gcdi_KPI(s_ed, s_uc), on = group_by_big)
  transform!(gcdi_KPI, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu)
  sort!(gcdi_KPI, :mu)
  return gcdi_KPI
end

function include_Δobjective_value(s_scalar, keys = [:day, :iteration], key = :objective_value, reference_configuration_key = :base_ramp_storage_envelopes_up_0_dn_0)
  # DEPRECATED
  ref_configuration = s_scalar[s_scalar.configuration .== reference_configuration_key, union([key], keys)]

  s_scalar = leftjoin(s_scalar, 
      rename(ref_configuration, key => :objective_value_ref_conf), 
      on = keys)

  s_scalar.Δobjective_value_ref_conf = s_scalar[!,key] .- s_scalar.objective_value_ref_conf
  s_scalar.Δobjective_value_relative_ref_conf = (s_scalar[!,key] .- s_scalar.objective_value_ref_conf)./s_scalar.objective_value_ref_conf
  return s_scalar
end

function calculate_adecuacy_gcd_KPI(gcdi_KPI, group_by = [:configuration, :day])
  gcd_KPI = outerjoin(
    combine(groupby(gcdi_KPI, group_by), [:LLD_h, :ENS_MWh] => ((x,y)->(LOLE = mean(x), EENS = mean(y))) => AsTable),  #TODO: change format
    combine(groupby(gcdi_KPI, group_by), [:CURD_h, :CUR_MWh] => ((x,y)->(CURE = mean(x), ECUR = mean(y))) => AsTable), #TODO: change format
    combine(groupby(gcdi_KPI, group_by), [:objective_value, :objective_value_uc, :OPEX, :OPEX_uc, :redispatch_cost, :LOL_cost, :LGEN_cost, :reserve_cost] .=> mean .=> [:EOV, :OV_uc, :EOPEX, :OPEX_uc, :E_redispatch_cost, :EENS_cost, :ELGEN_cost, :E_reserve_cost]), #:objective_value_uc changed to :OV_uc
    on=[:configuration, :day])
  if :LGEN_MWh in propertynames(gcdi_KPI) # ED
    gcd_KPI = outerjoin(gcd_KPI, combine(groupby(gcdi_KPI, group_by), :LGEN_MWh => mean => :ELGEN),on = [:configuration, :day])
  end
  transform!(gcd_KPI, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu)
  sort!(gcd_KPI, :mu)
  return gcd_KPI
end

function calculate_reserve_KPI(s_ed, s_uc, thres =.001)
  group_by = [:configuration, :day]
  group_by_big = [:configuration, :day, :iteration]
  KPI_reserve = leftjoin(
    s_ed.demand[!,union(group_by_big,[:hour, :demand_MW, :LOL_MW, :LGEN_MW])], 
    rename(s_uc.demand[!, union(group_by,[:hour, :demand_MW])], :demand_MW => :demand_uc_MW), 
    on = union(group_by,[:hour]))
  # leftjoin!(KPI_reserve, # used for calculation of Δenergy_reserve_up_MWh and Δenergy_reserve_dn_MWh
  #   combine(groupby(s_uc.storage, union(group_by,[:hour])), [:envelope_up_MWh, :envelope_down_MWh] .=> sum, renamecols = false),
  #   on=union(group_by,[:hour]))
  
  KPI_reserve.required_r_MW = (KPI_reserve.demand_MW .+ KPI_reserve.LOL_MW .- KPI_reserve.demand_uc_MW)
  KPI_reserve.required_r_relative = KPI_reserve.required_r_MW ./ KPI_reserve.demand_uc_MW

  KPI_reserve.redispatch_MW = KPI_reserve.LOL_MW - KPI_reserve.LGEN_MW  #- required_reserve.curtailment_MW
  KPI_reserve.redispatch_relative = KPI_reserve.redispatch_MW ./ KPI_reserve.demand_uc_MW
  KPI_reserve.redispatch_needed = (KPI_reserve.redispatch_MW.>=thres).*(KPI_reserve.required_r_MW.>=thres) .|| (KPI_reserve.redispatch_MW.<=-thres).*(KPI_reserve.required_r_MW.<=-thres) #not redispatch properly speaking but market deviation
  
  # TODO: revise the following ...
  KPI_reserve.required_r_up_MW =   KPI_reserve.required_r_MW.*(KPI_reserve.required_r_MW.>=0)
  KPI_reserve.required_r_dn_MW = -   KPI_reserve.required_r_MW.*(KPI_reserve.required_r_MW.<0)
  
  KPI_reserve = outerjoin(
    KPI_reserve,
    combine(groupby(s_ed.generation, union(group_by_big, [:hour])), :production_MW => sum => :production_MW),
    # combine(groupby(s_ed.storage, union(group_by_big, [:hour])), [:envelope_up_MWh, :envelope_down_MWh, :charge_MW, :discharge_MW, :SOE_MWh] .=> sum, renamecols = false),
    combine(groupby(s_ed.storage, union(group_by_big, [:hour])), [:charge_MW, :discharge_MW, :SOE_MWh] .=> sum, renamecols = false),
    on = union(group_by,[:iteration, :hour]))

  delivered = KPI_reserve.production_MW .+ KPI_reserve.discharge_MW .- KPI_reserve.charge_MW .- KPI_reserve.demand_uc_MW
  KPI_reserve.delivered_r_up_MW = delivered.*(delivered.>=0)
  KPI_reserve.delivered_r_dn_MW = -(delivered .+ KPI_reserve.LOL_MW).*(delivered.<0).*(KPI_reserve.required_r_MW.<0)
  KPI_reserve.delivered_r_up_ratio  = KPI_reserve.delivered_r_up_MW./KPI_reserve.required_r_up_MW
  KPI_reserve.delivered_r_dn_ratio  = KPI_reserve.delivered_r_dn_MW./KPI_reserve.required_r_dn_MW
  # end revise....
  
  # Calculation of Δenergy_reserve_up_MWh and Δenergy_reserve_dn_MWh
  # KPI_reserve.Δenergy_reserve_up_MWh = KPI_reserve.SOE_MWh - KPI_reserve.envelope_down_MWh
  # KPI_reserve.Δenergy_reserve_dn_MWh = KPI_reserve.envelope_up_MWh - KPI_reserve.SOE_MWh
  select!(KPI_reserve, Not([:charge_MW, :discharge_MW])) # cleaning Dataframe
  return sort(transform(KPI_reserve, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu), :mu)
end

function calculate_reserve_gcdi_KPI(KPI_reserve)
  group_by_big = [:configuration, :day, :iteration]
  gcdi_KPI_reserve  = combine(groupby(KPI_reserve, group_by_big),
    [:required_r_up_MW, :required_r_dn_MW, :delivered_r_up_MW, :delivered_r_dn_MW] .=> sum .=> [:required_r_up_MWh, :required_r_dn_MWh, :delivered_r_up_MWh, :delivered_r_dn_MWh])
  gcdi_KPI_reserve.delivered_r_up_ratio  = gcdi_KPI_reserve.delivered_r_up_MWh ./gcdi_KPI_reserve.required_r_up_MWh
  gcdi_KPI_reserve.delivered_r_dn_ratio  = gcdi_KPI_reserve.delivered_r_dn_MWh ./gcdi_KPI_reserve.required_r_dn_MWh
  leftjoin!(
    gcdi_KPI_reserve,
    # combine(groupby(KPI_reserve, group_by_big), [:Δenergy_reserve_up_MWh, :Δenergy_reserve_dn_MWh, :SOE_MWh].=> (x -> x[end]) .=> [:final_Δenergy_reserve_up_MWh, :final_Δenergy_reserve_dn_MWh, :final_SOE_MWh]),
    combine(groupby(KPI_reserve, group_by_big), [:SOE_MWh].=> (x -> x[end]) .=> [:final_SOE_MWh]),
    on = group_by_big)
  return sort(transform(gcdi_KPI_reserve, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu), :mu)
end

function calculate_reserve_gcd_KPI(gcdi_KPI_reserve)
  group_by = [:configuration, :day]
  gcd_KPI_reserve  = combine(groupby(gcdi_KPI_reserve, group_by), 
    # [:required_r_up_MWh, :required_r_dn_MWh, :delivered_r_up_MWh, :delivered_r_dn_MWh, :final_Δenergy_reserve_up_MWh, :final_Δenergy_reserve_dn_MWh, :final_SOE_MWh] .=> mean .=> [:E_required_r_up_MWh, :E_required_r_dn_MWh, :E_delivered_r_up_MWh, :E_deliverered_dn_MWh, :E_final_Δenergy_reserve_up_MWh, :E_final_Δenergy_reserve_dn_MWh, :E_final_SOE_MWh]
    [:required_r_up_MWh, :required_r_dn_MWh, :delivered_r_up_MWh, :delivered_r_dn_MWh, :final_SOE_MWh] .=> mean .=> [:E_required_r_up_MWh, :E_required_r_dn_MWh, :E_delivered_r_up_MWh, :E_deliverered_dn_MWh,  :E_final_SOE_MWh]
    )
  gcd_KPI_reserve.delivered_r_up_ratio  = gcd_KPI_reserve.E_delivered_r_up_MWh ./gcd_KPI_reserve.E_required_r_up_MWh
  gcd_KPI_reserve.delivered_r_dn_ratio  = gcd_KPI_reserve.E_deliverered_dn_MWh ./gcd_KPI_reserve.E_required_r_dn_MWh
  return sort(transform(gcd_KPI_reserve, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu), :mu)
end

function calculate_objective_function_gcdi_KPI(s_ed, s_uc)
  #TODO: add s_uc OPEX
  group_by = [:configuration, :day, :iteration]
  out =   combine(groupby(s_ed[:objective_function], group_by), 
    [:production_cost, :fixed_cost, :start_cost] => ((x,y,z) -> sum(skipmissing(x)) +  sum(skipmissing(y))) => :OPEX,
    [:LOL_cost, :LGEN_cost, :reserve_cost] .=> (x ->sum(skipmissing(x))) .=> [:LOL_cost, :LGEN_cost, :reserve_cost]
    )
  leftjoin!(out,
    combine(groupby(s_uc[:objective_function], [:configuration, :day]), 
      [:production_cost, :fixed_cost, :start_cost] => ((x,y,z) -> sum(skipmissing(x)) +  sum(skipmissing(y))) => :OPEX_uc,
      [:reserve_cost] .=> (x ->sum(skipmissing(x))) .=> [:reserve_cost_uc]),
      on = [:configuration, :day], 
  )
  out.objective_value = out.OPEX .+ out.LOL_cost .+ out.LGEN_cost .+ out.reserve_cost
  out.objective_value_uc = out.OPEX_uc .+ out.reserve_cost_uc
  out.redispatch_cost = out.OPEX .- out.OPEX_uc
  return sort(transform(out, :configuration .=> ByRow(x -> parse_configuration_to_mu(x)) .=> :mu), :mu)
end

function calculate_objective_function_gcd_KPI(gcdi_objective_funtion_KPI)
  #TODO: add s_uc OPEX
  group_by = [:configuration, :day, :mu]
  return  combine(groupby(gcdi_objective_funtion_KPI, group_by),
  [:OPEX, :redispatch_cost, :LOL_cost, :LGEN_cost, :reserve_cost] .=> mean .=> [:EOPEX, :E_redispatch_cost, :EENS_cost, :ELGEN_cost, :E_reserve_cost],
  )
end
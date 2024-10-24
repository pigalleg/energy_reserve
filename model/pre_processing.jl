using DataFrames
using CSV
using Distributions
using LinearAlgebra

G_DEFAULT_LOCATION = "./input/base_case"
G_NET_GENERAION_FULL_ID = "net_generation"
G_UC_DATA = "uc"



# --- start pre_processing ---

function variance(σ, ρ)
  # For autocorrelated errors X[t+1] = ρ*X[t] + (1-ρ^2)^(1/2)*N(0,σ[t])
  # then var(X[t]) = ρ^2*var(X[t-1]) + (1-ρ^2)*σ[t]^2, with var(X[0]) = σ[0]^2
  var = Vector{Float64}(undef,  length(σ))
  var[1] = σ[1]^2
  for i in 2:length(var)
    var[i] = ρ^2*var[i-1] + (1-ρ^2)* σ[i]^2
  end
  return var
end

function get_var(timeseries, ρ, p, margin = 0.1)
  σ = abs.(timeseries*margin./quantile(Normal(),p))
  σ = reshape(σ,(24,:))
  return mapslices(x->variance(x, ρ), σ, dims=1)
end

function to_GMT(df) # deprecated
  # Convert from GMT to GMT-8
  df.hour = mod.(df.hour .- 9, 8760) .+ 1
  sort!(df, :hour)
end

function generate_input_data(simulation_day, input_location = G_DEFAULT_LOCATION)
  gen_info, fuels, loads_df, gen_variable_info, storage_info = read_data(input_location)
  # gen_info, fuels, loads_df, gen_variable_info, storage_info = read_data()
  gen_df = pre_process_generators_data(gen_info, fuels)
  gen_df, loads_df, gen_variable_df  = pre_process_load_gen_variable(gen_df, loads_df, pre_process_gen_variable(gen_df, gen_variable_info))
  storage_df = pre_process_storage_data(storage_info)
  random_loads_df = read_random_demand(input_location)
  if !isnothing(simulation_day)
    return gen_df, 
      filter_day(simulation_day, loads_df),
      filter_day(simulation_day, gen_variable_df),
      storage_df,
      filter_day(simulation_day, random_loads_df)
  else
    return gen_df, loads_df, gen_variable_df, storage_df, random_loads_df
  end
end

function filter_periods(simulation_day, df)
  # deprecated
  T_period = (simulation_day*24+1):((simulation_day+1)*24)
  # Filtering data with timeseries according to T_period
  return df[in.(df.hour,Ref(T_period)),:]
end

function filter_day(simulation_day, df)
  if :day in propertynames(df)
    return df[in.(df.day,simulation_day),:]
  else
    T_period = ((simulation_day-1)*24+1):(((simulation_day-1)+1)*24)
    return df[in.(df.hour,Ref(T_period)),:]
  end
end


function filter_demand(expected_load, loads_to_filter, required_reserve)
  select = :day in propertynames(loads_to_filter) ? [:hour,:day] : [:hour]
  return transform(loads_to_filter, Not(select) .=> (x -> clamp.(x, expected_load.demand .- required_reserve.reserve_down_MW, expected_load.demand .+ required_reserve.reserve_up_MW)) .=> Not(select))
end

function read_data(input_location = G_DEFAULT_LOCATION; shift_timezone = false)
  input_uc_data_location = joinpath(input_location, G_UC_DATA)
  gen_info = CSV.read(joinpath(input_uc_data_location,"Generators_data.csv"), DataFrame)
  fuels = CSV.read(joinpath(input_uc_data_location,"Fuels_data.csv"), DataFrame)
  loads = CSV.read(joinpath(input_uc_data_location,"Demand.csv"), DataFrame)
  gen_variable = CSV.read(joinpath(input_uc_data_location,"Generators_variability.csv"), DataFrame)
  storage_info = CSV.read(joinpath(input_uc_data_location,"Storage_data.csv"), DataFrame)

  # rename all columns to lowercase (by convention)
  for f in [gen_info, fuels, loads, gen_variable, storage_info]
      rename!(f,lowercase.(names(f)))
  end
  if shift_timezone
    to_GMT(gen_variable)
    to_GMT(loads)
  end
  return gen_info, fuels, loads, identity.(gen_variable), storage_info
end

function pre_process_generators_data(gen_info,  fuels)
  # Keep columns relevant to our UC model
  select!(gen_info, 1:28) # columns 1:26
  gen_df = outerjoin(gen_info,  fuels, on = :fuel) # load in fuel costs and add to data frame
  rename!(gen_df, :cost_per_mmbtu => :fuel_cost)   # rename column for fuel cost
  gen_df.fuel_cost[ismissing.(gen_df[:,:fuel_cost])] .= 0

  # create "is_variable" column to indicate if this is a variable generation source (e.g. wind, solar):
  gen_df[!, :is_variable] .= false
  gen_df[in(["onshore_wind_turbine","small_hydroelectric","solar_photovoltaic", "net_generation"]).(gen_df.resource),
      :is_variable] .= true;

  # create full name of generator (including geographic location and cluster number)
  #  for use with variable generation dataframe
  gen_df.full_id = lowercase.(gen_df.region .* "_" .* gen_df.resource .* "_" .* string.(gen_df.cluster) .* ".0");

  # remove generators with no capacity (e.g. new build options that we'd use if this was capacity expansion problem)
  gen_df = gen_df[gen_df.existing_cap_mw .> 0,:]

  # net generation = -net_load for net_load < 0
  push!(gen_df, last(gen_df))
  last_ = nrow(gen_df)
  for k in names(gen_df)
    gen_df[last_, k] = ifelse(gen_df[last_, k] isa AbstractString, "", 0.0)
  end
  gen_df[last_, :r_id] = maximum(gen_df.r_id) + 1
  gen_df[last_, :resource] = G_NET_GENERAION_FULL_ID
  gen_df[last_, :full_id] = G_NET_GENERAION_FULL_ID
  gen_df[last_, :existing_cap_mw] = 0
  gen_df[last_, :var_om_cost_per_mwh] = 0
  gen_df[last_, :is_variable] = true
  gen_df[last_, :is_variable] = true
  gen_df[last_, :ramp_up_percentage] = 1
  gen_df[last_, :ramp_dn_percentage] = 1
  return identity.(gen_df)
end

function pre_process_storage_data(storage_info)
  df = copy(storage_info)
  df.full_id = lowercase.(storage_info.region .* "_" .* storage_info.resource .* "_" .* string.(storage_info.cluster) .* ".0");
  return df
end

function pre_process_load_gen_variable(gen_df, loads_df, gen_variable)
  # tranfers negative demand to generation
  filter = loads_df.demand.<0
  if !any(filter)
    net_generation = copy(loads_df)
    net_generation.generation .= 0
    installed_capacity = 0
    net_generation.cf .= 0
  else
    net_generation = loads_df[filter,:]
    net_generation.generation = - net_generation.demand
    loads_df[filter,:demand].=0
    installed_capacity = maximum(net_generation.generation)
    net_generation.cf = net_generation.generation./installed_capacity
  end
  net_generation.full_id .= G_NET_GENERAION_FULL_ID

  gen_df = copy(gen_df)
  gen_df[gen_df.full_id .== G_NET_GENERAION_FULL_ID, :existing_cap_mw] .= installed_capacity # Assumes that the element is already in the df

  gen_variable[gen_variable[!, :full_id] .== G_NET_GENERAION_FULL_ID,:cf] .= 0 # values reset to zero for re-iterations on the ED
  gen_variable[gen_variable[!, :full_id] .== G_NET_GENERAION_FULL_ID,:existing_cap_mw] .= installed_capacity # values reset to zero for re-iterations on the ED

  gen_variable = leftjoin(gen_variable, select(net_generation, Not([:demand,:generation])), on = [:hour, :full_id], makeunique = true)
  gen_variable = select(gen_variable, [:hour, :full_id, :r_id, :existing_cap_mw], [:cf_1,:cf] =>ByRow(coalesce) => [:cf])
  # gen_variable[gen_variable.full_id.== G_NET_GENERAION_FULL_ID, :existing_cap_mw] .= installed_capacity
  return gen_df, loads_df, sort(gen_variable,[:r_id,:hour])
end

function pre_process_gen_variable(gen_df, gen_variable_info)
  gen_variable_info[!,G_NET_GENERAION_FULL_ID] .= 0 # net generation = -net_load for net_load < 0
  aux = stack(gen_variable_info, Not(:hour), variable_name=:full_id, value_name=:cf)
  return innerjoin(aux,
    gen_df[gen_df.is_variable .== 1,[:r_id, :full_id, :existing_cap_mw]],
    on = :full_id)
end

function read_random_demand(input_location = G_DEFAULT_LOCATION)
  return CSV.read(joinpath(input_location, "ed", "random_demand.csv"), DataFrame)
end

function read_demand_scenarios(input_location)
  return CSV.read(joinpath(input_location, G_UC_DATA, "scenarios", "scenarios_demand.csv"), DataFrame)
end

function read_probability_scenarios(input_location)
  return CSV.read(joinpath(input_location, G_UC_DATA, "scenarios", "scenarios_probability.csv"), DataFrame)
end

function generate_scenarios_data(simulation_day, input_location = G_DEFAULT_LOCATION)
  return (demand = filter_day(simulation_day, read_demand_scenarios(input_location)), probability = read_probability_scenarios(input_location))
end

function read_reserve(input_location = G_DEFAULT_LOCATION)
  return CSV.read(joinpath(input_location, G_UC_DATA, "Reserve.csv"), DataFrame)
end


function read_parquet_and_convert(file)
  columns_to_symbol = [:configuration, :iteration]
  println(file)
  out = DataFrame(Parquet2.Dataset(file); copycols=false)
  for k in intersect(columns_to_symbol, propertynames(out))
    out[!, k] = Symbol.(out[!, k])
  end
  return out
end

# --- end pre_processing ---

function generate_configuration(key, storage_df; reserve=nothing, energy_reserve=nothing)
  function generate_envelope_configuration(μ_up, μ_dn, storage_df)
    return Dict(
        :ramp_constraints => true,
        :storage => storage_df,
        # :reserve => required_reserve,
        # :enriched_solution => true,
        :storage_envelopes => true,
        :μ_up => μ_up,
        :μ_dn => μ_dn)
  end
  envelope_config_key = match(r"base_ramp_storage_envelopes_up_(\w+)_dn_(\w+)", string(key))
  μ_up = parse(Float64, replace(envelope_config_key[1], "_" => "."))
  μ_dn = parse(Float64, replace(envelope_config_key[2], "_" => "."))
  out =  generate_envelope_configuration(μ_up, μ_dn, storage_df)
  if !isnothing(energy_reserve)
    out[:energy_reserve] = energy_reserve
  else 
    out[:reserve] = reserve
  end
  return out
end

function generate_reserves_old(loads, gen_variable, margin_percentage, baseload = 0)
  filter = gen_variable[!,:full_id] .== G_NET_GENERAION_FULL_ID
  net_gen = gen_variable[filter,:cf] .* gen_variable[filter,:existing_cap_mw]
  required_reserve = DataFrame(
    hour = loads[!,:hour],
    reserve_up_MW = baseload .+ (loads[!,:demand] .+ net_gen).*margin_percentage,
    reserve_down_MW = (loads[!,:demand] .+ net_gen).*margin_percentage)
  required_reserve = (required_reserve.>0).*required_reserve .- (required_reserve.<0).*required_reserve # negative to positive values
  return required_reserve
end

function generate_reserves(loads, gen_variable, ε, ρ; margin=0.1)
  filter = gen_variable[!,:full_id] .== G_NET_GENERAION_FULL_ID
  net_gen = gen_variable[filter,:cf] .* gen_variable[filter,:existing_cap_mw]
  p = 1-ε # 0.975
  μ = (loads[!,:demand] .+ net_gen)
  var =  get_var(μ , ρ, p, margin)
  σ = vec(sqrt.(var))
  normals = Normal.(0, σ)
  return DataFrame(
    hour = loads[!,:hour],
    reserve_up_MW = quantile.(normals, p),
    reserve_down_MW = -cquantile.(normals, p))
end

function generate_energy_reserves(loads, gen_variable, ε, ρ; margin=0.1)
  # Assumes that X[t] t = [1...24] are independent N(0,σ[t])
  # therefore Z[i,t] = sum_{τ=i}^t X[τ] is N(0,σ_Z[i,t]) with σ_Z[i,t] = (sum_{τ=i}^t σ[τ]^2)^1/2
  filter = gen_variable[!,:full_id] .== G_NET_GENERAION_FULL_ID
  net_gen = gen_variable[filter,:cf] .* gen_variable[filter,:existing_cap_mw]
  p = 1-ε
  # σ = quantile(Normal(),p)
  μ = (loads[!,:demand] .+ net_gen)
  var =  get_var(μ , ρ, p, margin) # calculates var[t] for each autocorrolelated X[t] 
  σ = zeros(size(var,1), size(var,1))
  t_H = diagm(loads[!,:hour])
  i_H = diagm(loads[!,:hour])
  for i in 1:size(σ,1), j in i:size(σ,2)
    σ[i,j] = sqrt(sum((var)[i:j])) # σ_Z[t] = (sum_{τ=i}^t σ[τ]^2)^1/2 with var[t] = σ[t]^2
    t_H[i,j] = t_H[i,i] + j-i
    i_H[i,j] = i_H[i,i]
  end
  normals = Normal.(0, σ)
  r_up = quantile.(normals, p)
  r_down = -cquantile.(normals, p)
  return DataFrame(
    i_hour = vcat([i_H[i,i:end] for i in 1:size(i_H,1)]...),
    t_hour = vcat([t_H[i,i:end] for i in 1:size(t_H,1)]...),
    reserve_up_MW =  vcat([r_up[i,i:end] for i in 1:size(r_up,1)]...),
    reserve_down_MW = vcat([r_down[i,i:end] for i in 1:size(r_down,1)]...)
  )
end


function generate_energy_reserves_deprecated(required_reserve)
  required_energy_reserve = [(row_1.hour, row_2.hour, row_1.reserve_up_MW*(row_1.hour == row_2.hour), row_1.reserve_down_MW*(row_1.hour == row_2.hour)) for row_1 in eachrow(required_reserve), row_2 in eachrow(required_reserve) if row_1.hour <= row_2.hour]
  required_energy_reserve = DataFrame(required_energy_reserve)
  required_energy_reserve = rename(required_energy_reserve, :1 => :i_hour, :2 => :t_hour, :3 => :reserve_up_MW, :4 => :reserve_down_MW,)
  return required_energy_reserve
end

function generate_energy_reserves_cumulative(required_reserve)
  required_energy_reserve_cumulated = [(row_1.hour, row_2.hour, sum(required_reserve[(required_reserve.hour .>= row_1.hour).&(required_reserve.hour .<= row_2.hour),:reserve_up_MW]), sum(required_reserve[(required_reserve.hour .>= row_1.hour).&(required_reserve.hour .<= row_2.hour),:reserve_down_MW])) for row_1 in eachrow(required_reserve), row_2 in eachrow(required_reserve) if row_1.hour <= row_2.hour]
  required_energy_reserve_cumulated = DataFrame(required_energy_reserve_cumulated)
  required_energy_reserve_cumulated = rename(required_energy_reserve_cumulated, :1 => :i_hour, :2 => :t_hour, :3 => :reserve_up_MW, :4 => :reserve_down_MW,)
  return required_energy_reserve_cumulated
end
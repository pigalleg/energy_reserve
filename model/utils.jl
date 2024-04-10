using DataFrames
using CSV
using Random
using Distributions

G_DEFAULT_LOCATION = "./input/base_case"
G_NET_GENERAION_FULL_ID = "net_generation"
function to_GMT(df)
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

  loads_multi_df, gen_variable_multi_df, random_loads_multi_df = filter_periods(simulation_day, loads_df, gen_variable_df, random_loads_df)
  return gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df
end

function filter_periods(simulation_day, loads_df, gen_variable_df, random_loads_df)
  T_period = (simulation_day*24+1):((simulation_day+1)*24)
  # Filtering data with timeseries according to T_period
  gen_variable_multi_df = gen_variable_df[in.(gen_variable_df.hour,Ref(T_period)),:]
  loads_multi_df = loads_df[in.(loads_df.hour,Ref(T_period)),:]
  random_loads_multi_df =  random_loads_df[in.(random_loads_df.hour,Ref(T_period)),:]

  return loads_multi_df, gen_variable_multi_df, random_loads_multi_df
end

function read_data(input_location = G_DEFAULT_LOCATION)
  input_data_location = joinpath(input_location, "uc_data")
  gen_info = CSV.read(joinpath(input_data_location,"Generators_data.csv"), DataFrame)
  fuels = CSV.read(joinpath(input_data_location,"Fuels_data.csv"), DataFrame)
  loads = CSV.read(joinpath(input_data_location,"Demand.csv"), DataFrame)
  gen_variable = CSV.read(joinpath(input_data_location,"Generators_variability.csv"), DataFrame)
  storage_info = CSV.read(joinpath(input_data_location,"Storage_data.csv"), DataFrame)

  # Rename all columns to lowercase (by convention)
  for f in [gen_info, fuels, loads, gen_variable, storage_info]
      rename!(f,lowercase.(names(f)))
  end
  to_GMT(gen_variable)
  to_GMT(loads)
  return gen_info, fuels, loads, identity.(gen_variable), storage_info
end



function pre_process_generators_data(gen_info,  fuels)
  # Keep columns relevant to our UC model 
  select!(gen_info, 1:26) # columns 1:26
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

# Functions to randomize demand. Not used by the rest of the code
function create_random_demand(loads, nb_samples, input_location = G_DEFAULT_LOCATION)
  # example: create_random_demand(loads, 1000)
  out = transform(loads, :demand => ByRow(x ->  [rand(Normal(x, abs(x)*0.1/2)) for i in 1:nb_samples]) => ["demand_$i" for i in 1:nb_samples])
  CSV.write(joinpath(input_location, "demand","random_demand.csv"), out)
end

function read_random_demand(input_location = G_DEFAULT_LOCATION)
  return CSV.read(joinpath(input_location, "demand", "random_demand.csv"), DataFrame)
end

function main()
  gen_info, fuels, loads, gen_variable_info, storage_info = read_data()
  create_random_demand(loads, 1000)
end

function generate_configurations(storage_df, required_reserve, required_energy_reserve, required_energy_reserve_cumulated)
  configs = (
      base = (
          ramp_constraints = false,
          # storage = storage_df,
          # reserve = required_reserve,
          # energy_reserve = required_energy_reserve,
          enriched_solution = true,
          # storage_envelopes = false
      ),
      base_ramp = (
          ramp_constraints = true,
          # storage = storage_df,
          # reserve = required_reserve,
          # energy_reserve = required_energy_reserve,
          enriched_solution = true,
          # storage_envelopes = false
      ),
      base_ramp_reserve = (
          ramp_constraints = true,
          # storage = storage_df,
          reserve = required_reserve,
          # energy_reserve = required_energy_reserve,
          enriched_solution = true,
          # storage_envelopes = false
      ),
      base_ramp_energy_reserve = (
          ramp_constraints = true,
          # storage = storage_df,
          # reserve = required_reserve,
          energy_reserve = required_energy_reserve,
          enriched_solution = true,
          # storage_envelopes = false
      ),
      base_ramp_storage = (
          ramp_constraints = true,
          storage = storage_df,
          # reserve = required_reserve,
          # energy_reserve = required_energy_reserve,
          enriched_solution = true,
          # storage_envelopes = false
      ),
      base_ramp_storage_reserve = (
          ramp_constraints = true,
          storage = storage_df,
          reserve = required_reserve,
          # energy_reserve = required_energy_reserve,
          enriched_solution = true,
          # storage_envelopes = false
      ),
      base_ramp_storage_envelopes = (
          ramp_constraints = true,
          storage = storage_df,
          reserve = required_reserve,
          # energy_reserve = required_energy_reserve,
          enriched_solution = true,
          storage_envelopes = true
      ),
      base_ramp_storage_energy_reserve = (
          ramp_constraints = true,
          storage = storage_df,
          # reserve = required_reserve,
          energy_reserve = required_energy_reserve,
          enriched_solution = true,
          storage_envelopes = false,
      ),
      base_ramp_storage_energy_reserve_cumulated = (
          ramp_constraints = true,
          storage = storage_df,
          # reserve = required_reserve,
          energy_reserve = required_energy_reserve_cumulated,
          enriched_solution = true,
          storage_envelopes = false,
      ),
  )
  return Dict(pairs(configs))
end

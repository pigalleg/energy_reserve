using DataFrames
using CSV
using Random
using Distributions

function to_GMT(df)
  # Convert from GMT to GMT-8
  df.hour = mod.(df.hour .- 9, 8760) .+ 1
  sort!(df, :hour)
end

function read_data(base_location = ".")
  # input_data_location = joinpath("input", "uc_data")
  input_data_location = joinpath(base_location,"input", "uc_data")
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

function pre_process_gen_variable(gen_df, gen_variable_info)
  aux = stack(gen_variable_info, Not(:hour), variable_name=:full_id, value_name=:cf)
  return innerjoin(aux,
    gen_df[gen_df.is_variable .== 1,[:r_id, :full_id, :existing_cap_mw]],
    on = :full_id)
end

function pre_process_generators_data(gen_info,  fuels)
  # Keep columns relevant to our UC model 
  select!(gen_info, 1:26) # columns 1:26
  gen_df = outerjoin(gen_info,  fuels, on = :fuel) # load in fuel costs and add to data frame
  rename!(gen_df, :cost_per_mmbtu => :fuel_cost)   # rename column for fuel cost
  gen_df.fuel_cost[ismissing.(gen_df[:,:fuel_cost])] .= 0

  # create "is_variable" column to indicate if this is a variable generation source (e.g. wind, solar):
  gen_df[!, :is_variable] .= false
  gen_df[in(["onshore_wind_turbine","small_hydroelectric","solar_photovoltaic"]).(gen_df.resource),
      :is_variable] .= true;

  # create full name of generator (including geographic location and cluster number)
  #  for use with variable generation dataframe
  gen_df.full_id = lowercase.(gen_df.region .* "_" .* gen_df.resource .* "_" .* string.(gen_df.cluster) .* ".0");

  # remove generators with no capacity (e.g. new build options that we'd use if this was capacity expansion problem)
  gen_df = gen_df[gen_df.existing_cap_mw .> 0,:]
  return identity.(gen_df)
end

function pre_process_storage_data(storage_info)
  df = copy(storage_info)
  df.full_id = lowercase.(storage_info.region .* "_" .* storage_info.resource .* "_" .* string.(storage_info.cluster) .* ".0");
  return df
end

# Functions to randomize demand. Not used by the rest of the code
function create_random_demand(loads, nb_samples, base_location = ".")
  out = transform(loads, :demand => ByRow(x ->  [rand(Normal(x, x*0.1/2)) for i in 1:nb_samples]) => ["demand_$i" for i in 1:nb_samples])
  CSV.write(joinpath(base_location, "input", "demand","random_demand.csv"), out)
end

function read_random_demand(base_location = ".")
  return CSV.read(joinpath(base_location, "input", "demand", "random_demand.csv"), DataFrame)
end

function main()
  gen_info, fuels, loads, gen_variable_info, storage_info = read_data()
  create_random_demand(loads, 1000)
end

function generate_configurations(required_energy_reserve, required_energy_reserve_cumulated)
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

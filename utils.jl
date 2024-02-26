using DataFrames
using CSV

function to_GMT(df)
  # Convert from GMT to GMT-8
  df.hour = mod.(df.hour .- 9, 8760) .+ 1
  sort!(df, :hour)
end

function read_data()
  datadir = joinpath("input", "uc_data")
  gen_info = CSV.read(joinpath(datadir,"Generators_data.csv"), DataFrame)
  fuels = CSV.read(joinpath(datadir,"Fuels_data.csv"), DataFrame)
  loads = CSV.read(joinpath(datadir,"Demand.csv"), DataFrame)
  gen_variable = CSV.read(joinpath(datadir,"Generators_variability.csv"), DataFrame)
  storage_info = CSV.read(joinpath(datadir,"Storage_data.csv"), DataFrame)

  # Rename all columns to lowercase (by convention)
  for f in [gen_info, fuels, loads, gen_variable, storage_info]
      rename!(f,lowercase.(names(f)))
  end
  to_GMT(gen_variable)
  to_GMT(loads)
  return gen_info, fuels, loads, gen_variable, storage_info
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
  return gen_df
end

function pre_process_storage_data(storage_info)
  df = copy(storage_info)
  df.full_id = lowercase.(storage_info.region .* "_" .* storage_info.resource .* "_" .* string.(storage_info.cluster) .* ".0");
  return df
end





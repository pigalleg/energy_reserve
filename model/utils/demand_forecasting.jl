using DataFrames
using CSV
using Random
using Distributions
include("../utils.jl")

G_INSAMPLE_ERRORS = "insample_errors.csv"

# Functions to randomize demand. Not used by the rest of the code
function create_random_demand(loads, nb_samples, output_file)
  # example: create_random_demand(loads, 1000)
  out = transform(loads, :demand => ByRow(x ->  [rand(Normal(x, abs(x)*0.1/2)) for i in 1:nb_samples]) => ["scenario_$i" for i in 1:nb_samples])
  select!(out, Not(:demand))
  CSV.write(output_file, out)
end

function create_correlated_random_demand(loads, nb_samples, input_location)
  in_sample_errors = CSV.read(G_INSAMPLE_ERRORS, DataFrame)
  # @infiltrate

end

function main()
  input_location = "../../input/base_case_increased_storage_energy_v3.1"
  output_file = "scenarios_demand.csv"
  gen_info, fuels, loads, gen_variable_info, storage_info = read_data(input_location)
  create_random_demand(loads, 1000, output_file)
end

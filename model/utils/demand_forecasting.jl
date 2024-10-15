using DataFrames
using CSV
using Random
using Distributions
using Infiltrator
include("../pre_processing.jl")

# Functions to randomize demand. Not used by the rest of the code
# function create_random_demand(loads, nb_samples, output_file)
#   margin = 0.1
#   q = 0.975
#   factor = quantile(Normal(),q)
#   out = transform(loads, :demand => ByRow(x ->  [rand(Normal(x, abs(x)*margin/factor)) for i in 1:nb_samples]) => ["scenario_$i" for i in 1:nb_samples])
#   select!(out, Not(:demand))
#   CSV.write(output_file, out)
# end


# main("../../input/ercot_brownfield_expansion_v1.0")
function main(input_location)
  # input_location = "../../input/ercot_brownfield_expansion_v1.0"
  gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(nothing, input_location)
  ρ = 0.8
  p = 0.975
  # create_autocorrelated_demand(loads_multi_df, ρ, p) #-->"../../input/base_case_increased_storage_energy_v4.1"
  create_reserve(random_loads_multi_df, p)
end


function test(loads, q=0.975)
  margin = 0.1
  # epsilon = 0.025
  factor = quantile(Normal(),q)
  # normales = Normal.(loads.demand, abs.(loads.demand*margin./factor))
  normales = Normal.(0, abs.(loads.demand*margin./factor))
  quantiles =  quantiles = quantile.(normales, q)
  p = cdf.(normales, quantiles)
  # output should be
  # q = 0.975 --> e = 0.025
  # q = 0.025 --> e = 0.025
end

function check()
  input_location = "../../input/base_case_increased_storage_energy_v4.1"
  gen_info, fuels, loads, gen_variable_info, storage_info = read_data(input_location)
  ρ = 0.8
  p = 0.975
  d = loads.demand[25:48]
  n  = 10000000
  var = get_var(d, ρ, p)
  σ = vec(sqrt.(var))
  normals = Normal.(0, σ)
  percentile_1 = quantile.(normals, p)
  errors = create_autocorrelated_errors(d, ρ, p; n_errors=n)
  percentile_2 = mapslices(x->quantile(x, p), errors, dims=2)
  @show (percentile_1.-percentile_2)./percentile_1
  @infiltrate
end


function check2()
  input_location = "../../input/base_case_increased_storage_energy_v4.1"
  gen_info, fuels, loads, gen_variable_info, storage_info = read_data(input_location)
  ρ = 0.8
  p = 0.975
  d = loads.demand[1:48]
  # d = ones(48) # if d (standard deviation) is constant for each day, then corr are the same for each days
  # d[25:46].=2
  n  = 100000
  errors = create_autocorrelated_errors(d, ρ, p; n_errors=n)
  # errors = reshape(errors,(24,2,n))
  errors = reshape(errors,(24,n,2))
  covar = get_covar(d, ρ, p)
  corr = stack([correlation(covar[:,:,j],ρ) for j in 1:size(covar,3)])
  mvnormal = mapslices(x ->MvNormal(zeros(24),x), covar, dims = [1,2])
  cov_ = cov.(mvnormal)
  cor_ = cor.(mvnormal)  
  @show cov_[:,:,2][1].-covar[:,:,2]
  @show cor_[:,:,2][1].-corr[:,:,2]
  empirical_cov = mapslices(x ->cov(x, dims = 2), errors, dims = [1,2])
  empirical_cor = mapslices(x ->cor(x, dims = 2), errors, dims = [1,2])
  # @show empirical_cor.-cor_
  @show (empirical_cov[:,:,2] .-covar[:,:,2])
  @show (empirical_cor[:,:,2] .-corr[:,:,2])
  @infiltrate
end

function check3()
  input_location = "../../input/base_case_increased_storage_energy_v4.1"
  gen_info, fuels, loads, gen_variable_info, storage_info = read_data(input_location)
  n = 10000
  ρ = 0.8
  p = 0.975
  errors = create_autocorrelated_errors(loads.demand, ρ, p; n_errors = n)
  errors = reshape(errors,(24,365,n))
  # covar = get_covar(d, ρ, p)
  # mvnormal = mapslices(x ->MvNormal(zeros(24),x), covar, dims = [1,2])
  # cov_ = cov.(mvnormal)[1]
  # cor_ = cor.(mvnormal)[1]  
  # empirical_cov = cov(errors, dims = 2)
  empirical_cor =  mapslices(x ->cor(x, dims = 2), errors, dims = [1,2])
  @infiltrate
end


function create_autocorrelated_demand(loads, ρ, p)
  output_file = "random_demand.csv"
  errors = create_autocorrelated_errors(loads.demand, ρ, p)
  demand =  DataFrame(errors.+ loads.demand, ["demand_$i" for i in 1:size(errors,2)])
  insertcols!(demand,1, :demand => loads.demand)
  insertcols!(demand,1, :hour => loads.hour)
  CSV.write(output_file, demand)
end

function create_autocorrelated_errors(d, ρ, p; n_errors = 10000)
  covar = get_covar(d, ρ, p)
  mvnormal = mapslices(x ->MvNormal(zeros(24),x), covar, dims = [1,2])
  errors = rand.(mvnormal, n_errors) # n_samples
  errors = vcat(errors...)
  return errors
end

function create_reserve(random_loads, p)
  select_ = :day in propertynames(random_loads) ? [:hour,:day] : [:hour]
  errors = transform(random_loads, Not(select_) .=> (x -> x.- random_loads.demand), renamecols = false)
  select!(errors, Not(:demand))
  errors = stack(errors, Not(select_))
  grouped = groupby(errors, select_)
  reserve = combine(grouped,
    :value => (x -> quantile(x,p)) => :reserve_up_MW,
    :value => (x -> -quantile(x,1-p)) => :reserve_down_MW,
    )
  CSV.write("Reserve.csv", reserve)
  @show reserve
end


function get_covar(timeseries, ρ, p)
  var = get_var(timeseries, ρ, p)
  covar = stack([covariance(y,ρ) for y in eachcol(var)])
  return covar
end

function covariance(var, ρ)
  # For autocorrelated errors X[t+1] = ρ*X[t] + (1-ρ^2)^(1/2)*N(0,σ[t])
  # then cov(X[t],X[t+1]) = ρ*cov(X[t],X[t])
  covar = Matrix(Diagonal(var))
  for i in 1:(size(covar, 1)),j = (i+1):(size(covar, 2))
    covar[i,j] =  ρ * covar[i,j-1]
    covar[j,i] = covar[i,j]
  end
  return covar
end

function correlation(cov, ρ)
  # by definition of corr
  cor = copy(cov)
  for i in 1:(size(cor, 1)),j = (i):(size(cor, 2))
    cor[i,j] =  cov[i,j]/sqrt(cov[i,i]*cov[j,j])
    cor[j,i] = cor[i,j]
  end 
  return cor
end

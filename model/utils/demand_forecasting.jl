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



function generate_autocorrelated_demand_file(input_location, ρ)
  println("Generating autocorrelated demand file with ρ = $ρ")
  # generate_autocorrelated_demand_file("../../input/ercot_brownfield_expansion_v1.0")
  loads_df = CSV.read(joinpath(input_location, "uc", "Demand.csv"), DataFrame)
  p = 0.975
  create_autocorrelated_demand(loads_df, ρ, p, input_location) #-->"../../input/base_case_increased_storage_energy_v4.1"
end

function generate_reserve_file(input_location)
  println("Generating reserve file")  
  random_load = CSV.read(joinpath(input_location, "ed", "random_demand.csv"), DataFrame)
  p = 0.975
  create_reserve(random_load, p, input_location)
end

function generate_energy_reserve_file(input_location)
  println("Generating energy reserve file")
  random_load = CSV.read(joinpath(input_location, "ed", "random_demand.csv"), DataFrame)
  p = 0.975
  create_energy_reserve(random_load, p, input_location)
end


function create_autocorrelated_demand(loads, ρ, p, output_location = nothing)
  errors = create_autocorrelated_errors(loads.demand, ρ, p)
  demand =  DataFrame(errors.+ loads.demand, ["demand_$i" for i in 1:size(errors,2)])
  insertcols!(demand,1, :demand => loads.demand)
  insertcols!(demand,1, :hour => loads.hour)
  if :day in propertynames(loads)
    insertcols!(demand,1, :day => loads.day)
  end

  if !isnothing(output_location)
    output_file = joinpath(output_location, "ed", "random_demand.csv")
  else  
    output_file = "random_demand.csv"
  end
  println("Saving random demand file at $output_file")
  CSV.write(output_file, demand)
end

function create_autocorrelated_errors(d, ρ, p; n_errors = 1000)
  covar = get_covar(d, ρ, p)
  mvnormal = mapslices(x ->MvNormal(zeros(24),x), covar, dims = [1,2])
  errors = rand.(mvnormal, n_errors) # n_samples
  errors = vcat(errors...)
  return errors
end

function create_reserve(random_loads, p, output_location = nothing) 
  select_ = :day in propertynames(random_loads) ? [:hour,:day] : [:hour]
  errors = transform(random_loads, Not(select_) .=> (x -> x.- random_loads.demand), renamecols = false)
  select!(errors, Not(:demand))
  errors = stack(errors, Not(select_))
  grouped = groupby(errors, select_)
  reserve = combine(grouped,
    :value => (x -> quantile(x,p)) => :reserve_up_MW,
    :value => (x -> -quantile(x,1-p)) => :reserve_down_MW,
    )
  if !isnothing(output_location)
    output_file = joinpath(output_location, "uc", "Reserve.csv")
  else  
    output_file = "Reserve.csv"
  end 
  println("Saving reserve file at $output_file")
  CSV.write(output_file, reserve)
end

function create_energy_reserve(random_loads, p, output_location = nothing)
  # It needs :day in random_loads
  # tuples_ij(hour) = [(i_hour = i, t_hour = t) for i in hour, t in hour if i<=t] 
  # percentile_ij_(hour, value) = [sum((i.<=hour .* hour.<=t).*value) for i in unique(hour), t in unique(hour) if i<=t] 
  function aux_(i_hour, df, p)
    df_ = []
    reserve_up_ = []
    reserve_down_ = []
    i_ = []
    t_= []
    for t_hour in df.hour if t_hour >= i_hour
          push!(i_, i_hour)
          push!(t_, t_hour)
        if t_hour == i_hour
          push!(df_, df[df.hour.== t_hour, Not([:hour,:day])])
        else
          push!(df_, df[df.hour.== t_hour, Not([:hour,:day])] .+ last(df_))
        end
        push!(reserve_up_, quantile(last(df_)[1,:], p))
        push!(reserve_down_, -quantile(last(df_)[1, :], 1-p))
      end
    end
    return (i_hour = i_, t_hour = t_, reserve_up_MW = reserve_up_, reserve_down_MW = reserve_down_)
  end
  errors = transform(random_loads, Not([:hour,:day]) .=> (x -> x.- random_loads.demand), renamecols = false)
  select!(errors, Not(:demand))
  energy_reserve = combine(groupby(errors, [:day]), AsTable(:) => (x -> [aux_(i, DataFrame(x), p) for i in x.hour]) => AsTable)
  energy_reserve = combine(groupby(energy_reserve, [:day]), Not(:day) .=> (x -> reduce(vcat,x)),  renamecols = false)
  if !isnothing(output_location)
    output_file = joinpath(output_location, "uc", "Energy reserve.csv")
  else  
    output_file = "Energy reserve.csv"
  end 
  println("Saving energy reserve file at $output_file")
  CSV.write(output_file, energy_reserve)
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

function check1()
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
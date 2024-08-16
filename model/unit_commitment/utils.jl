function remove_variable_constraint(model, key, delete_ = true)
  # Applies for constraints and variables
  println("Removing $key...")
  if !haskey(model, key)
      println("variable not in model")
      return
  end
  if delete_ delete.(model, model[key]) end # Constraints must be deleted also
  unregister(model, key)
end

# TODO change gen_variable => gen_varialbe_df, loads => loads_df
function create_generators_sets(gen_df)
  # Thermal resources for which unit commitment constraints apply
  G_thermal = gen_df[gen_df[!,:up_time] .> 0,:r_id] 
      
  # Non-thermal resources for which unit commitment constraints do NOT apply 
  G_nonthermal = gen_df[gen_df[!,:up_time] .== 0,:r_id]
  
  # Variable renewable resources
  G_var = gen_df[gen_df[!,:is_variable] .== 1,:r_id]
  
  # Non-variable (dispatchable) resources
  G_nonvar = gen_df[gen_df[!,:is_variable] .== 0,:r_id]
  
  # Non-variable and non-thermal resources
  G_nt_nonvar = intersect(G_nonvar, G_nonthermal)
  # Note that G_nt_var = G_var

  # Set of all generators (above are all subsets of this)~
  G = gen_df.r_id

  return G, G_thermal, G_nonthermal, G_var, G_nonvar, G_nt_nonvar
end

function create_time_sets(loads)
  return loads.hour, loads.hour[1:end-1]
end

function get_sets(gen_df, scenarios)
  G, G_thermal, G_nonthermal, G_var, G_nonvar, G_nt_nonvar = create_generators_sets(gen_df)
  T, T_red = create_time_sets(scenarios.demand)
  Σ = create_scenarios_sets(scenarios.probability)
  return sets = (G = G,
    G_thermal = G_thermal,
    G_nonthermal = G_nonthermal,
    G_var = G_var,
    G_nonvar = G_nonvar,
    G_nt_nonvar = G_nt_nonvar,
    T = T,
    T_red = T_red,
    Σ = Σ)
end

function create_storage_sets(storage)
  return storage.r_id
end

function create_scenarios_sets(scenarios_probability)
    return scenarios_probability.scenario
end

function convert_to_indexed_vector(value, T)
  if value isa Number
      return Dict(T .=> fill(value, length(T)))
  else
      return Dict(T .=> value)
  end
end
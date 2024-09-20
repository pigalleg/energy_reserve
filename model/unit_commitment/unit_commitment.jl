using DataFrames
include("./deterministic.jl")
include("./stochastic.jl")
# include("../utils.jl")
# include("../post_processing.jl")

function construct_deterministic_unit_commitment(gen_df, loads, gen_variable, mip_gap, storage, ramp_constraints; kwargs...)
    println("Constructing DUC...")
    reserve = get(kwargs, :reserve, nothing)
    energy_reserve = get(kwargs, :energy_reserve, nothing)
    storage_envelopes = get(kwargs, :storage_envelopes, false)
    storage_link_constraint =  get(kwargs, :storage_link_constraint, false)
    μ_up = get(kwargs, :μ_up, 1)
    μ_dn = get(kwargs, :μ_dn, 1)
    VRESERVE = get(kwargs, :value_reserve, 1e-6)
    bidirectional_storage_reserve = get(kwargs, :bidirectional_storage_reserve, true)
    thermal_reserve = get(kwargs, :thermal_reserve, false)
    naive_envelopes = get(kwargs, :naive_envelopes, false)
    sets =  get_sets(gen_df, loads)
    uc = DUC(gen_df, loads, gen_variable, mip_gap)
    if !isnothing(storage)
        println("Adding storage...")
        add_storage(uc, storage, loads, gen_df, sets)
    end
    if ramp_constraints
        println("Adding ramp constraints...")   
        add_ramp_constraints(uc, gen_df, sets)
    end
    if !isnothing(reserve) 
        println("Adding reserve constraints...")
        add_reserve_constraints(uc, reserve, loads, gen_df, storage, bidirectional_storage_reserve, storage_envelopes, naive_envelopes, thermal_reserve, μ_up, μ_dn, VRESERVE, sets)
    end
    if !isnothing(energy_reserve)
        println("Adding energy reserve constraints...")
        add_energy_reserve_constraints(uc, energy_reserve, loads, gen_df, storage, storage_link_constraint, thermal_reserve, VRESERVE, sets)
    end
    return uc
end

function construct_stochastic_unit_commitment(gen_df, gen_variable, mip_gap, storage, ramp_constraints, scenarios, expected_min_SOE, VLOL, VLGEN)
    println("Constructing SUC...")
    uc = SUC(gen_df, gen_variable, scenarios, mip_gap, VLOL, VLGEN)
    if !isnothing(storage)
        println("Adding storage...")
        add_storage_s(uc, storage, scenarios, get_sets(gen_df, scenarios.demand, scenarios.probability), expected_min_SOE) # TODO: This function is meant to be used within SUC
    end
    if ramp_constraints
        println("Adding ramp constraints...")   
        add_ramp_constraints_s(uc, gen_df, get_sets(gen_df, scenarios.demand, scenarios.probability)) # TODO: This function is meant to be used within SUC
    end
    return uc
end


function construct_unit_commitment(gen_df, loads, gen_variable, scenarios; kwargs...)
    storage = get(kwargs, :storage, nothing)
    ramp_constraints = get(kwargs, :ramp_constraints, false)
    mip_gap = get(kwargs, :mip_gap, 1e-8)
    stochastic = get(kwargs, :stochastic, false)
    if !stochastic
        return construct_deterministic_unit_commitment(gen_df, loads, gen_variable, mip_gap, storage, ramp_constraints; kwargs...)
    else
        return construct_stochastic_unit_commitment(gen_df, gen_variable, mip_gap, storage, ramp_constraints, scenarios, get(kwargs, :expected_min_SOE, false), get(kwargs, :VLOL, 1e4), get(kwargs, :VLGEN, 0))
    end
end

function solve_unit_commitment(gen_df, loads, gen_variable, scenarios = nothing; kwargs...)
    reference_solution =  get(kwargs, :reference_solution, nothing)
    uc = construct_unit_commitment(gen_df, loads, gen_variable, scenarios; kwargs...)
    # relax_integrality(uc)
    # include("./debugging_ignore.jl")
    # set_optimizer_attribute(model, "OutputFlag", 1)
    if !isnothing(reference_solution)
        uc = generate_alternative_model(uc, reference_solution)
    end
    # save_model_to_file(uc,"uc")
    optimize!(uc)
    if !is_solved_and_feasible(uc)
        @infiltrate
        # list = get_conflicting_constraints(uc)
    end
    return uc
end
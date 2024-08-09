using DataFrames
include("./deterministic.jl")
include("./stochastic.jl")
include("../utils.jl")
include("../post_processing.jl")

function construct_deterministic_unit_commitment(gen_df, loads, gen_variable, mip_gap, storage, ramp_constraints; kwargs...)
    reserve = get(kwargs, :reserve, nothing)
    energy_reserve = get(kwargs, :energy_reserve, nothing)
    storage_envelopes = get(kwargs, :storage_envelopes, false)
    storage_link_constraint =  get(kwargs, :storage_link_constraint, false)
    μ_up = get(kwargs, :μ_up, 1)
    μ_dn = get(kwargs, :μ_dn, 1)
    VRESERVE = get(kwargs, :value_reserve, 1e-6)
    bidirectional_storage_reserve = get(kwargs, :bidirectional_storage_reserve, true)
    thermal_reserve = get(kwargs, :thermal_reserve, false)
    uc = unit_commitment(gen_df, loads, gen_variable, mip_gap)
    if !isnothing(storage)
        println("Adding storage...")
        add_storage(uc, storage, loads, gen_df)
    end
    if ramp_constraints
        println("Adding ramp constraints...")   
        add_ramp_constraints(uc, gen_df)
    end
    if !isnothing(reserve) 
        println("Adding reserve constraints...")
        add_reserve_constraints(uc, reserve, loads, gen_df, storage, bidirectional_storage_reserve, storage_envelopes, thermal_reserve, μ_up, μ_dn, VRESERVE)
    end
    if !isnothing(energy_reserve)
        println("Adding energy reserve constraints...")
        add_energy_reserve_constraints(uc, energy_reserve, loads, gen_df, storage, storage_link_constraint)
    end
    return uc
end

function construct_stochastic_unit_commitment(gen_df, loads, gen_variable, mip_gap, storage, ramp_constraints, scenarios)
    uc = unit_commitment(gen_df, loads, gen_variable, scenarios, mip_gap)
    if !isnothing(storage)
        println("Adding storage...")
        add_storage(uc, storage, scenarios)
    end
    if ramp_constraints
        println("Adding ramp constraints...")   
        add_ramp_constraints(uc, gen_df, scenarios)
    end
    return uc
end


function construct_unit_commitment(gen_df, loads, gen_variable, scenarios; kwargs...)
    storage = get(kwargs, :storage, nothing)
    ramp_constraints = get(kwargs, :ramp_constraints, false)
    mip_gap = get(kwargs, :mip_gap, 1e-8)
    stochastic = get(kwargs, :stochastic, false)
    println("Constructing UC...")
    if !stochastic
        return construct_deterministic_unit_commitment(gen_df, loads, gen_variable, mip_gap, storage, ramp_constraints; kwargs...)
        
    else
        return construct_stochastic_unit_commitment(gen_df, loads, gen_variable, mip_gap, storage, ramp_constraints, scenarios)
    end
end

function solve_unit_commitment(gen_df, loads, gen_variable, scenarios; kwargs...)
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
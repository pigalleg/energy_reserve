using Infiltrator
include("./model/pre_processing.jl")
include("./model/post_processing.jl")
include("./model/unit_commitment/unit_commitment.jl")
include("./model/economic_dispatch.jl")
# include("./notebooks/plotting.jl")
include("./notebooks/processing.jl")

# __revise_mode__ = :eval
# ENV["COLUMNS"]=120 # Set so all columns of DataFrames and Matrices are displayed
function plot_results(solution)
    supply, demand = calculate_supply_demand(solution)
    p1 = plot_fieldx_by_fieldy(supply, :production_MW, :resource)
    p2 = plot_fieldx_by_fieldy(demand, :demand_MW, :resource)
    p3 = plot()
    p4 = plot()
    p5 = plot()
    p6 = plot()
    has_reserve = false
    if haskey(solution,:energy_reserve)
        solution_reserve = solution.energy_reserve[solution.energy_reserve.hour.==solution.energy_reserve.hour_i,:]
        has_reserve = true
    elseif haskey(solution,:reserve)
        solution_reserve = copy(solution.reserve)
        has_reserve = true
    end
    if has_reserve
        reserve = calculate_reserve(solution_reserve, required_reserve)
        p3 = plot_reserve_by_fieldy(reserve, :reserve_up_MW, :resource)
        p4 = plot_reserve_by_fieldy(reserve, :reserve_down_MW, :resource)
    end

    if haskey(solution,:storage) & has_reserve
        battery_reserve = calculate_battery_reserve(solution.storage, solution_reserve)
        p5 = plot_battery_reserve_(battery_reserve, :reserve_up_MW_eff)
        p6 = plot_battery_reserve_(battery_reserve, :reserve_down_MW_eff)
    end
    [
        p1 p2
        p3 p4
        p5 p6
    ]
end

function generate_multipliers_configurations(μs)
    mu_to_string(x) = isinteger(x) ? string(Int(x)) : replace(string(x), "." => "_")
    return [Symbol("base_ramp_storage_envelopes_up_$(mu_to_string(μ))_dn_$(mu_to_string(μ))") for μ in μs]
end

G_day = 7 # 68
G_RESERVE = 0.1
G_input_folder = "./input/base_case_increased_storage_energy_v3.1"
# G_REMOVE_RESERVE_CONSTRAINTS = true
# G_CONSTRAIN_DISPATCH = true
# G_MAX_ITERATIONS = 100
# G_VRESERVE = 1e-6
# G_REMOVE_VARIABLES_FROM_OBJECTIVE = false


config = (
    ramp_constraints = true,
    # energy_reserve = required_energy_reserve,
    # energy_reserve = required_energy_reserve_cumulated,
    enriched_solution = true,
    
    # μ_up = 1,
    # μ_dn = 1,
)

# function load_configuration(storage_df, required_reserve)
#     config = (
#         mip_gap = G_MIP_GAP,
#         ramp_constraints = true,
#         storage = storage_df,
#         # reserve = required_reserve,
#         # energy_reserve = required_energy_reserve,
#         # energy_reserve = required_energy_reserve_cumulated,
#         enriched_solution = true,
#         # storage_envelopes = true,
#         # μ_up = 1,
#         # μ_dn = 1,
#     )
#     return config
# end

function load_deterministic_data(day, input_folder, reserve = G_RESERVE)
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(day, input_folder)
    required_reserve = generate_reserves(loads_multi_df, gen_variable_multi_df, reserve)
    random_loads_multi_df = filter_demand(loads_multi_df, random_loads_multi_df, required_reserve)
    return gen_df, loads_multi_df, random_loads_multi_df, gen_variable_multi_df, storage_df, required_reserve
end

function load_scenarios(day, input_folder, loads_multi_df, required_reserve)
    scenarios = generate_scenarios_data(day, input_folder)
    scenarios = (demand = filter_demand(loads_multi_df, scenarios.demand, required_reserve), probability = scenarios.probability)
    return scenarios
end

function duc(;kwargs...)
    input_folder = get(kwargs, :input_folder, G_input_folder)
    day = get(kwargs, :day, G_day)
    gen_df, loads_multi_df, random_loads_multi_df, gen_variable_multi_df, storage_df, required_reserve = load_deterministic_data(day, input_folder)
    return solve_unit_commitment(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df;
        storage = storage_df,
        # reserve = required_reserve,
        # storage_envelopes = true,
        energy_reserve = generate_energy_reserves(required_reserve),
        storage_link_constraint = true,        
        config...
        )
end

function suc(;kwargs...)
    input_folder = get(kwargs, :input_folder, G_input_folder)
    day = get(kwargs, :day, G_day)
    expected_min_SOE = get(kwargs, :expected_min_SOE, false)
    gen_df, loads_multi_df, random_loads_multi_df, gen_variable_multi_df, storage_df, required_reserve = load_deterministic_data(day, input_folder)
    scenarios = load_scenarios(day, input_folder, loads_multi_df, required_reserve)
    return solve_unit_commitment(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df,
        scenarios;
        storage = storage_df,
        stochastic = true,
        expected_min_SOE = expected_min_SOE,
        config...
        )

end

function ed()
    solution  = solve_economic_dispatch(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df;
        config...
        )
    return solution
end

function ed_multi_demand(iteration_key)
    solution  = solve_economic_dispatch(
        gen_df,
        random_loads_multi_df,
        gen_variable_multi_df;
        config...
        )
    return NamedTuple(k => filter(:iteration => ==(iteration_key), v) for (k,v) in zip(keys(solution), solution))
end

function uc_net_demand()
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_day, "./input/net_demand_case")
    required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, G_RESERVE)
    solution  = solve_unit_commitment(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df;
        config...
        )
    return solution
end

function ed_multi_demand_net_demand()
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_day, "./input/net_demand_case")
    required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, G_RESERVE)
    solution  = solve_economic_dispatch(
        gen_df,
        random_loads_multi_df,
        gen_variable_multi_df;
        config...
        )
    return solution
end

function generate_ed_solutions_(days, configurations; kwargs...)

    function get_reference_configuration(k, configurations)
        i = first(findall(x->x == k , configurations))
        return i > 1 ?  getindex(configurations, i-1) : nothing
    end

    input_folder = get(kwargs, :input_folder, G_input_folder)
    output_folder = get(kwargs, :output_folder, "./output")
    write = get(kwargs, :write, true)
    reserve = get(kwargs, :reserve, 0.1)
    alternative_solution = get(kwargs, :alternative_solution, false)

    add_config = Dict(
        # :max_iterations => get(kwargs, :max_iterations, 100),
        # :constrain_dispatch => get(kwargs, :constrain_dispatch, true),
        # :value_reserve => get(kwargs, :value_reserve, 1e-6),
        # :remove_variables_from_objective => get(kwargs, :remove_variables_from_objective, false),
        # :VLOL => get(kwargs, :VLOL, 1e4),
        # :VLGEN => get(kwargs, :VLGEN, 0),
        # :thermal_reserve =>  get(kwargs, :thermal_reserve, false),
        # :bidirectional_storage_reserve => get(kwargs, :bidirectional_storage_reserve, true),
        :constrain_SOE_by_envelopes => get(kwargs, :constrain_SOE_by_envelopes, false),
        # :constrain_dispatch_by_multipliers => get(kwargs, :constrain_dispatch_by_multipliers, false),
        # :naive_envelopes => get(kwargs, :naive_envelopes, false),
        # :variables_to_constrain => get(kwargs, :variables_to_constrain, [GEN, CH, DIS])
        # :stochastic => get(kwargs, :stochastic, false),
    )
    # configurations = vcat(configurations, [:base_ramp_storage_energy_reserve_cumulated])
    s_uc = Dict()
    s_ed = Dict()
    for day in days, k in configurations

        gen_df, loads_multi_df, random_loads_multi_df, gen_variable_multi_df, storage_df, required_reserve = load_deterministic_data(day, input_folder, reserve)
        config = merge(add_config, generate_configuration(k, storage_df, required_reserve)[k])
        
        k_reference = get_reference_configuration(k, configurations)
        if !isnothing(k_reference) & alternative_solution # if reference_solution is added, both uc and ed are will be solved with alternative model
            config = merge((reference_solution = s_uc[(day,k_reference)],), config)
        end
        uc = solve_unit_commitment(
            gen_df,
            loads_multi_df,
            gen_variable_multi_df;
            config...
        )
        
        s_uc[(day,k)] = get_model_solution(uc, gen_df, loads_multi_df, gen_variable_multi_df; config...)
        s_ed[(day,k)] = solve_economic_dispatch_get_solution(
            uc,
            gen_df,
            random_loads_multi_df,
            gen_variable_multi_df;
            config...
        )
    end
    s_ed = merge_solutions(s_ed, [:day, :configuration])
    s_uc = merge_solutions(s_uc, [:day, :configuration])

    s_uc = Dict(pairs(s_uc))
    if haskey(s_uc, :energy_reserve) s_uc[:reserve] = vcat(s_uc[:reserve], s_uc[:energy_reserve][s_uc[:energy_reserve].hour.==s_uc[:energy_reserve].hour_i,:][:,Not(:hour_i)]) end
    s_uc = NamedTuple(s_uc)

    if write
        if !isdir(output_folder) mkdir(output_folder) end
        folder_path = joinpath(output_folder,"n_$(join(days,"-"))")
        solution_to_parquet(s_uc, "s_uc", folder_path)
        solution_to_parquet(s_ed, "s_ed", folder_path)
    end
    return s_uc, s_ed
end

function generate_ed_solutions(;days, μs, kwargs...)
    for day in days
        generate_ed_solutions_([day], generate_multipliers_configurations(μs); kwargs...)
    end
    generate_post_processing_KPI_files(get(kwargs, :output_folder, nothing))
end

function merge_ed_solutions(solution_folders, folder_path, read = true, write = false)
    if read
        keys = [:demand, :generation, :storage, :reserve, :energy_reserve, :scalar]
        s_uc = [parquet_to_solution("s_uc", joinpath(folder_path, s)) for s in solution_folders]
        s_ed = [parquet_to_solution("s_ed", joinpath(folder_path, s)) for s in solution_folders]
        s_uc = NamedTuple(k => vcat([s[k] for s in s_uc if haskey(s, k)]...) for k in keys)
        s_ed = NamedTuple(k => vcat([s[k] for s in s_ed if haskey(s, k)]...) for k in keys)
    end
    if write
        name = "n_$(replace(join(solution_folders, "-"), "n_" =>""))"
        solution_to_parquet(s_uc, "s_uc", joinpath(folder_path, name))
        solution_to_parquet(s_ed, "s_ed", joinpath(folder_path, name))
    end
    return s_uc, s_ed
end

function generate_post_processing_KPI_files(folder_path, folders_to_read_ = nothing, save = true)
    function chunk_list_custom(arr, chunk_size)
        return [arr[i:min(i + chunk_size - 1, end)] for i in 1:chunk_size:length(arr)]
    end

    function KPI_df_list(solution_folders, folder_path)
        println("Calculating adecuacy KPIS...")
        s_uc, s_ed = merge_ed_solutions(solution_folders, folder_path)
        gcdi_KPI_adequacy = calculate_adecuacy_gcdi_KPI(s_ed, s_uc)
        gcd_KPI_adequacy = calculate_adecuacy_gcd_KPI(gcdi_KPI_adequacy)
        KPI_reserve = calculate_reserve_KPI(s_ed, s_uc)
        gcdi_KPI_reserve = calculate_reserve_gcdi_KPI(KPI_reserve)
        gcd_KPI_reserve = calculate_reserve_gcd_KPI(gcdi_KPI_reserve)
        println("...done)")
        return [gcdi_KPI_adequacy, gcd_KPI_adequacy, KPI_reserve, gcdi_KPI_reserve, gcd_KPI_reserve]

    end
    chunk_size = 1
    values = [[],[],[],[],[]]

    folders_to_read = last.(splitpath.(filter(isdir, readdir(folder_path; join = true))))
    out_name = "all"
    if !isnothing(folders_to_read_)
        folders_to_read = intersect(folders_to_read_, folders_to_read)
        out_name = join(folders_to_read, "_")
    end
    keys = [:gcdi_KPI_adequacy, :gcd_KPI_adequacy, :KPI_reserve, :gcdi_KPI_reserve, :gcd_KPI_reserve]
    values = vcat.([KPI_df_list(folders, folder_path) for folders in chunk_list_custom(folders_to_read, chunk_size)]...)
    out = NamedTuple(k => v for (k,v) in zip(keys, values))
   
    if save
        solution_to_parquet(out, out_name, folder_path)
    end
    return out
end
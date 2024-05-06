include("./model/utils.jl")
include("./model/economic_dispatch.jl")
include("./notebooks/plotting.jl")
include("./notebooks/processing.jl")
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

G_N = 7 # 68
G_RESERVE = 0.1
G_MIP_GAP = 10^(-8)
G_REMOVE_RESERVE_CONSTRAINTS = true
G_CONSTRAIN_DISPATCH = true
G_MAX_ITERATIONS = 10

gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_N)
# random_loads_multi_df = random_loads_multi_df[!, [:hour, :demand, :demand_21]]
# random_loads_multi_df.demand_21 .= random_loads_multi_df.demand_21*0.5
required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, G_RESERVE)

config = (
    ramp_constraints = true,
    storage = storage_df,
    reserve = required_reserve,
    # energy_reserve = required_energy_reserve,
    # energy_reserve = required_energy_reserve_cumulated,
    enriched_solution = true,
    storage_envelopes = true,
    μ_up = 1,
    μ_dn = 1,
)
ed_config = (
    remove_reserve_constraints = G_REMOVE_RESERVE_CONSTRAINTS,
    max_iterations = G_MAX_ITERATIONS,
    constrain_dispatch = G_CONSTRAIN_DISPATCH
)
config = merge(config, ed_config)


function uc()
    solution  = solve_unit_commitment(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df,
        G_MIP_GAP;
        config...
        )
    return solution
end

function ed()
    solution  = solve_economic_dispatch(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df,
        G_MIP_GAP;
        config...
        )
    return solution
end

function ed_multi_demand(iteration_key)
    solution  = solve_economic_dispatch(
        gen_df,
        random_loads_multi_df,
        gen_variable_multi_df,
        G_MIP_GAP;
        config...
        )
    return NamedTuple(k => filter(:iteration => ==(iteration_key), v) for (k,v) in zip(keys(solution), solution))
end

function uc_net_demand()
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_N, "./input/net_demand_case")
    required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, G_RESERVE)
    solution  = solve_unit_commitment(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df,
        G_MIP_GAP;
        config...
        )
    return solution
end

function ed_multi_demand_net_demand()
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_N, "./input/net_demand_case")
    required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, G_RESERVE)
    solution  = solve_economic_dispatch(
        gen_df,
        random_loads_multi_df,
        gen_variable_multi_df,
        G_MIP_GAP;
        config...
        )
    return solution
end

function generate_ed_solutions_(days, configurations; kwargs...)
    
    # if isnothing(configurations)
    #     μs = [(0,0), (0.25,0.25), (0.5,0.5), (0.75, 0.75), (1, 1)]
    #     configurations = [Symbol("base_ramp_storage_envelopes_up_$(replace(string(μ_up), "." => "_"))_dn_$(replace(string(μ_dn), "." => "_"))") for (μ_up, μ_dn) in μs]
        
    # end
    input_folder = get(kwargs, :input_folder, "./input/base_case")
    output_folder = get(kwargs, :output_folder, "./output")
    max_iterations = get(kwargs, :max_iterations, 100)
    write = get(kwargs, :write, true)
    reserve = get(kwargs, :reserve, G_RESERVE)
    constrain_dispatch = get(kwargs, :constrain_dispatch, G_CONSTRAIN_DISPATCH)
    remove_reserve_constraints = get(kwargs, :remove_reserve_constraints, G_REMOVE_RESERVE_CONSTRAINTS)

    ed_config = Dict(:max_iterations => max_iterations, :constrain_dispatch => constrain_dispatch, :remove_reserve_constraints => remove_reserve_constraints)
    configurations = vcat(configurations, [:base_ramp_storage_energy_reserve_cumulated])
    s_uc = Dict()
    s_ed = Dict()
    
    for day in days, k in configurations
        gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(day, input_folder)
        required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, reserve)
        config = Dict(k => merge(v, ed_config) for (k,v) in generate_configuration(k, storage_df, required_reserve, required_energy_reserve, required_energy_reserve_cumulated))
        s_uc[(day,k)] = solve_unit_commitment(
            gen_df,
            loads_multi_df,
            gen_variable_multi_df,
            G_MIP_GAP;
            config[k]...
        )
        s_ed[(day,k)] = solve_economic_dispatch(
            gen_df,
            random_loads_multi_df,
            gen_variable_multi_df,
            G_MIP_GAP;
            config[k]...
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
        println("Calculating adecuacy KPIS...)")
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
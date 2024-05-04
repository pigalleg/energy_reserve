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

G_N = 7 # 68
G_RESERVE = 0.1
G_MIP_GAP = 0.00000001
G_REMOVE_RESERVE_CONSTRAINTS = true
G_CONSTRAIN_DISPATCH = true
G_MAX_ITERATIONS = 10

gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_N)
random_loads_multi_df = random_loads_multi_df[!, [:hour, :demand, :demand_21]]
random_loads_multi_df.demand_21 .= random_loads_multi_df.demand_21*0.5
# println(names(random_loads_multi_df))
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

function ed_multi_demand()
    solution  = solve_economic_dispatch(
        gen_df,
        random_loads_multi_df,
        gen_variable_multi_df,
        G_MIP_GAP;
        config...
        )
    return solution
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

function generate_ed_solutions_(days, input_folder; max_iterations = 100, configurations = nothing, output_folder = ".", write = true, reserve = G_RESERVE)
    ed_config = Dict(:max_iterations => max_iterations, :constrain_dispatch => G_CONSTRAIN_DISPATCH, :remove_reserve_constraints => G_REMOVE_RESERVE_CONSTRAINTS)
    if isnothing(configurations)
        μs = [(0,0), (0.25,0.25), (0.5,0.5), (0.75, 0.75), (1, 1)]
        configurations = [Symbol("base_ramp_storage_envelopes_up_$(replace(string(μ_up), "." => "_"))_dn_$(replace(string(μ_dn), "." => "_"))") for (μ_up, μ_dn) in μs]
        other_configs = [:base_ramp_storage_energy_reserve_cumulated]
        configurations = vcat(configurations, other_configs)
    end
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

function generate_ed_solutions(days, input_folder; max_iterations = 100, configurations = nothing, output_folder = ".", write = true, reserve = G_RESERVE)
    # generate_ed_solutions([15, 45, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350], "./input/base_case", output_folder = "./output/solutions_v3")
    for day in days
        generate_ed_solutions_([day], input_folder,  max_iterations = max_iterations, configurations = configurations,  output_folder = output_folder, write = write, reserve = reserve)
    end
end

function merge_ed_solutions(solution_folders, folder_path)
    # merge_ed_solutions(["n_15", "n_45", "n_75", "n_106", "n_136", "n_167", "n_197", "n_228"], joinpath(".","output", "solutions_v1.2"))
    read = true
    write = true
    if read
        keys = [:demand, :generation, :storage, :reserve, :energy_reserve, :scalar]
        s_uc = [parquet_to_solution("s_uc", joinpath(folder_path, s)) for s in solution_folders]
        s_ed = [parquet_to_solution("s_ed", joinpath(folder_path, s)) for s in solution_folders]
        s_uc = NamedTuple(k => vcat([s[k] for s in s_uc]..., cols = :union) for k in keys)
        s_ed = NamedTuple(k => vcat([s[k] for s in s_ed]..., cols = :union) for k in keys)
    end
    # @infiltrate
    if write
        name = "n_$(replace(join(solution_folders, "-"), "n_" =>""))"
        solution_to_parquet(s_uc, "s_uc", joinpath(folder_path, name))
        solution_to_parquet(s_ed, "s_ed", joinpath(folder_path, name))
    end
    return s_uc, s_ed
end
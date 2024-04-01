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

n=100
gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(n)
required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, 0.05, 0)

config = (
    ramp_constraints = true,
    storage = storage_df,
    # reserve = required_reserve,
    energy_reserve = required_energy_reserve,
    # energy_reserve = required_energy_reserve_cumulated,
    enriched_solution = true,
    storage_envelopes = true
)
ed_config = (
    remove_reserve_constraints = true,
    max_iterations = 10,
)
config = merge(config, ed_config)


function main_uc()
    solution  = solve_unit_commitment(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df,
        0.0001;
        config...
        )
    plot_results(solution)
end

function main_ed()
    solution  = solve_economic_dispatch(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df,
        0.0001;
        config...
        )
        plot_results(solution)
end

function main_ed_multi_demand()
    solution  = solve_economic_dispatch(
        gen_df,
        random_loads_multi_df,
        gen_variable_multi_df,
        0.0001;
        config...
        )
        plot_results(solution)
end

function main_uc_net_demand()
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(n, "./input/net_demand_case")
    required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, 0.05, 0)
    solution  = solve_unit_commitment(
        gen_df,
        loads_multi_df,
        gen_variable_multi_df,
        0.0001;
        config...
        )
    plot_results(solution)
end

# main_uc()
# main_ed()
# main_ed_multi_demand()
main_uc_net_demand()

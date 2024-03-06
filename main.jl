# using Plots
# using PrettyTables
# using VegaLite
# using Debugger
# using Infiltrator
include("./utils.jl")
include("./unit_commitment.jl")
include("./plotting.jl")
# ENV["COLUMNS"]=120 # Set so all columns of DataFrames and Matrices are displayed

gen_info, fuels, loads, gen_variable_info, storage_info = read_data()
gen_df = pre_process_generators_data(gen_info, fuels)
storage_df = pre_process_storage_data(storage_info)

# A spring day
n=100
T_period = (n*24+1):((n+1)*24)

# High solar case: 3,500 MW
gen_variable = pre_process_gen_variable(gen_df, gen_variable_info)

# No thermal generation
# gen_df[gen_df.resource .== "natural_gas_fired_combustion_turbine",:existing_cap_mw] .= 0.1
# gen_df[gen_df.resource .== "natural_gas_fired_combined_cycle",:existing_cap_mw] .= 0.1

# Filtering data with timeseries according to T_period
gen_variable_multi = gen_variable[in.(gen_variable.hour,Ref(T_period)),:];
loads_multi = loads[in.(loads.hour,Ref(T_period)),:];

required_reserve = DataFrame(
    hour = loads[in.(loads.hour, Ref(T_period)), :hour],
    reserve_up_MW = 300 .+ loads[in.(loads.hour,Ref(T_period)), :demand].*0.05,
    reserve_down_MW = loads[in.(loads.hour, Ref(T_period)), :demand].*0.05)

required_energy_reserve = [(row_1.hour, row_2.hour, row_1.reserve_up_MW*(row_1.hour == row_2.hour), row_1.reserve_down_MW*(row_1.hour == row_2.hour)) for row_1 in eachrow(required_reserve), row_2 in eachrow(required_reserve) if row_1.hour <= row_2.hour]
required_energy_reserve = DataFrame(required_energy_reserve)
required_energy_reserve = rename(required_energy_reserve, :1 => :i_hour, :2 => :t_hour, :3 => :reserve_up_MW, :4 => :reserve_down_MW,)


required_energy_reserve_cumulated = [(row_1.hour, row_2.hour, sum(required_reserve[(required_reserve.hour .>= row_1.hour).&(required_reserve.hour .<= row_2.hour),:reserve_up_MW]), sum(required_reserve[(required_reserve.hour .>= row_1.hour).&(required_reserve.hour .<= row_2.hour),:reserve_down_MW])) for row_1 in eachrow(required_reserve), row_2 in eachrow(required_reserve) if row_1.hour <= row_2.hour]
required_energy_reserve_cumulated = DataFrame(required_energy_reserve_cumulated)
required_energy_reserve_cumulated = rename(required_energy_reserve_cumulated, :1 => :i_hour, :2 => :t_hour, :3 => :reserve_up_MW, :4 => :reserve_down_MW,)
;

function main()
    solution  = solve_unit_commitment(
        gen_df,
        loads_multi,
        gen_variable_multi,
        0.01,
        ramp_constraints = true,
        storage = storage_df,
        # reserve = required_reserve,
        energy_reserve = required_energy_reserve_cumulated,
        enriched_solution = true,
        # storage_envelopes = true
        )
    supply, demand = calculate_supply_demand(solution)
    if haskey(solution,:energy_reserve)
        solution_reserve = solution.energy_reserve[solution.energy_reserve.hour.==solution.energy_reserve.hour_i,:]
    else
        solution_reserve = copy(solution.reserve)
    end
    reserve = calculate_reserve(solution_reserve, required_reserve)
    # @infiltrate
    battery_reserve = calculate_battery_reserve(solution.storage, solution_reserve)
    # @infiltrate
    [
        plot_fieldx_by_fieldy(supply, :production_MW, :resource) plot_fieldx_by_fieldy(demand, :demand_MW, :resource)
        plot_reserve_by_fieldy(reserve, :reserve_up_MW, :resource) plot_reserve_by_fieldy(reserve, :reserve_down_MW, :resource)
        # plot_fieldx_by_fieldy(solution_reserve[solution_reserve.resource .== "battery",:], :reserve_up_MW, :r_id)
        plot_fieldx_by_fieldy(battery_reserve, :envelope_up_MW_eff, :resource) plot_fieldx_by_fieldy(demand, :demand_MW, :resource)
    ]
end
main()
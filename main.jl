# using Plots
# using PrettyTables
# using VegaLite
using Debugger
using Infiltrator
include("./utils.jl")
include("./unit_commitment.jl")
# ENV["COLUMNS"]=120 # Set so all columns of DataFrames and Matrices are displayed

gen_info, fuels, loads, gen_variable_info, storage_info = read_data()
gen_df = pre_process_generators_data(gen_info, fuels)
storage_df = pre_process_storage_data(storage_info)

# A spring day
n=100
T_period = (n*24+1):((n+1)*24)

# High solar case: 3,500 MW
gen_df_sens = copy(gen_df)
gen_df_sens[gen_df_sens.resource .== "solar_photovoltaic",
    :existing_cap_mw] .= 3500
gen_variable = pre_process_gen_variable(gen_df_sens, gen_variable_info)

# Filtering data with timeseries according to T_period
gen_variable_multi = gen_variable[in.(gen_variable.hour,Ref(T_period)),:];
loads_multi = loads[in.(loads.hour,Ref(T_period)),:];

requested_reserve = DataFrame(
    hour = loads[in.(loads.hour, Ref(T_period)), :hour],
    reserve_up_MW = 300 .+ loads[in.(loads.hour,Ref(T_period)), :demand].*0.05,
    reserve_down_MW = loads[in.(loads.hour, Ref(T_period)), :demand].*0.05)

requested_energy_reserve = [(row_1.hour, row_2.hour, row_1.reserve_up_MW*(row_1.hour == row_2.hour), row_1.reserve_down_MW*(row_1.hour == row_2.hour)) for row_1 in eachrow(requested_reserve), row_2 in eachrow(requested_reserve) if row_1.hour <= row_2.hour]
requested_energy_reserve = DataFrame(requested_energy_reserve)
requested_energy_reserve = rename(requested_energy_reserve, :1 => :i_hour, :2 => :t_hour, :3 => :reserve_up_MW, :4 => :reserve_down_MW,)

function main()
    solution = solve_unit_commitment(
        gen_df_sens,
        loads_multi,
        gen_variable_multi,
        0.001,
        ramp_constraints = true,
        storage = storage_df,
        # reserve = requested_reserve,
        energy_reserve = requested_energy_reserve,
        
        enriched_solution = true)
    # @infiltrate
    println(solution.generation)
    println(solution.storage)
    println(solution.demand)
    # println(solution.re)
    # 
end

main()
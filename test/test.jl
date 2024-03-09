

include("../utils.jl")
include("../unit_commitment.jl")
include("../economic_dispatch.jl")
# include("../plotting.jl")
using CSV
OBJECTIVE_VALUE = :objective_value

reference =  CSV.read(joinpath("./test","reference.csv"), DataFrame)

gen_info, fuels, loads, gen_variable_info, storage_info = read_data()
gen_df = pre_process_generators_data(gen_info, fuels)
storage_df = pre_process_storage_data(storage_info)

# A spring day
n=100
T_period = (n*24+1):((n+1)*24)

# High solar case: 3,500 MW
gen_variable = pre_process_gen_variable(gen_df, gen_variable_info)

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

configs = (
    base = (
        ramp_constraints = false,
        # storage = storage_df,
        # reserve = required_reserve,
        # energy_reserve = required_energy_reserve,
        enriched_solution = true,
        # storage_envelopes = false
    ),
    base_ramp = (
        ramp_constraints = true,
        # storage = storage_df,
        # reserve = required_reserve,
        # energy_reserve = required_energy_reserve,
        enriched_solution = true,
        # storage_envelopes = false
    ),
    base_ramp_reserve = (
        ramp_constraints = true,
        # storage = storage_df,
        reserve = required_reserve,
        # energy_reserve = required_energy_reserve,
        enriched_solution = true,
        # storage_envelopes = false
    ),
    base_ramp_energy_reserve = (
        ramp_constraints = true,
        # storage = storage_df,
        # reserve = required_reserve,
        energy_reserve = required_energy_reserve,
        enriched_solution = true,
        # storage_envelopes = false
    ),
    base_ramp_storage = (
        ramp_constraints = true,
        storage = storage_df,
        # reserve = required_reserve,
        # energy_reserve = required_energy_reserve,
        enriched_solution = true,
        # storage_envelopes = false
    ),
    base_ramp_storage_reserve = (
        ramp_constraints = true,
        storage = storage_df,
        reserve = required_reserve,
        # energy_reserve = required_energy_reserve,
        enriched_solution = true,
        # storage_envelopes = false
    ),
    base_ramp_storage_envelopes = (
        ramp_constraints = true,
        storage = storage_df,
        reserve = required_reserve,
        # energy_reserve = required_energy_reserve,
        enriched_solution = true,
        storage_envelopes = true
    ),
    base_ramp_storage_energy_reserve = (
        ramp_constraints = true,
        storage = storage_df,
        # reserve = required_reserve,
        energy_reserve = required_energy_reserve,
        enriched_solution = true,
        storage_envelopes = false,
    ),
    base_ramp_storage_energy_reserve_cumulated = (
        ramp_constraints = true,
        storage = storage_df,
        # reserve = required_reserve,
        energy_reserve = required_energy_reserve_cumulated,
        enriched_solution = true,
        storage_envelopes = false,
    ),
)

function test1()
    out = DataFrame(Dict(string(k) => solve_unit_commitment(
            gen_df,
            loads_multi,
            gen_variable_multi,
            0.001;
            v...).objective_value for (k,v) in zip(keys(configs), configs)
    ))
    # @infiltrate
    println(out == reference)
end

function test2()
    mip_gap = 0.0001
    # configs_ = (base = configs[:base],)
    out = [(solve_unit_commitment(
            gen_df,
            loads_multi,
            gen_variable_multi,
            mip_gap;
            v...)[OBJECTIVE_VALUE],
        solve_economic_dispatch_single_demand(
            gen_df,
            loads_multi,
            gen_variable_multi,
            mip_gap;
            v...)[OBJECTIVE_VALUE]) for (k,v) in zip(keys(configs), configs)
    ]
    out = DataFrame(Config = collect(keys(configs)), UC= getindex.(out,1), EC = getindex.(out,2))
    out[!,:delta_percentage] .= (out.UC .- out.EC)./out.UC*100
    out[!,:delta_percentage_loq_mip_gap] .= out.delta_percentage .<=mip_gap
    println(out)
end
# test1()
test2()
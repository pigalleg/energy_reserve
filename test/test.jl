

include("../utils.jl")
include("../unit_commitment.jl")
include("../economic_dispatch.jl")
include("../processing.jl")
# include("../plotting.jl")
using CSV
OBJECTIVE_VALUE = :objective_value
SCALAR = :scalar

reference =  CSV.read(joinpath("./test","reference.csv"), DataFrame)

gen_info, fuels, loads, gen_variable_info, storage_info = read_data()
gen_df = pre_process_generators_data(gen_info, fuels)
gen_variable = pre_process_gen_variable(gen_df, gen_variable_info)
storage_df = pre_process_storage_data(storage_info)
random_loads_df = read_random_demand()

# A spring day
n=100
T_period = (n*24+1):((n+1)*24)

# Filtering data with timeseries according to T_period
gen_variable_multi = gen_variable[in.(gen_variable.hour,Ref(T_period)),:];
loads_multi = loads[in.(loads.hour,Ref(T_period)),:]
random_loads_multi =  random_loads_df[in.(random_loads_df.hour,Ref(T_period)),:];

required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi, 0.05, 300)
configs = generate_configurations(required_energy_reserve, required_energy_reserve_cumulated)

ed_config = (
    remove_reserve_constraints = false,
)

function test1()
    out = DataFrame(Dict(string(k) => solve_unit_commitment(
            gen_df,
            loads_multi,
            gen_variable_multi,
            0.001;
            v...).scalar[1, OBJECTIVE_VALUE] for (k,v) in configs
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
            v...).scalar[1, OBJECTIVE_VALUE],
            solve_economic_dispatch(
            gen_df,
            loads_multi,
            gen_variable_multi,
            mip_gap;
            merge(v, ed_config)...).scalar[1, OBJECTIVE_VALUE] ) for (k,v) in configs
    ]
    out = DataFrame(Config = collect(keys(configs)), UC= getindex.(out,1), EC = getindex.(out,2))
    out[!,:delta_percentage] .= (out.UC .- out.EC)./out.UC*100
    out[!,:delta_percentage_loq_mip_gap] .= out.delta_percentage .<=mip_gap
    println(out)
end
# test1()
# test2()
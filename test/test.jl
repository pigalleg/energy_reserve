

include("../model/utils.jl")
include("../model/unit_commitment.jl")
include("../model/economic_dispatch.jl")
include("../notebooks/processing.jl")
# include("../plotting.jl")
using CSV
OBJECTIVE_VALUE = :objective_value
SCALAR = :scalar


# gen_info, fuels, loads, gen_variable_info, storage_info = read_data("./input/net_demand_case")
# # gen_info, fuels, loads, gen_variable_info, storage_info = read_data()
# gen_df = pre_process_generators_data(gen_info, fuels)
# loads, gen_variable, gen_df = pre_process_load_gen_variable(loads, gen_df, gen_variable_info)
# storage_df = pre_process_storage_data(storage_info)
# random_loads_df = read_random_demand()

# # A spring day

# T_period = (n*24+1):((n+1)*24)

# # Filtering data with timeseries according to T_period
# gen_variable_multi = gen_variable[in.(gen_variable.hour,Ref(T_period)),:]
# loads_multi = loads[in.(loads.hour,Ref(T_period)),:]
# random_loads_multi =  random_loads_df[in.(random_loads_df.hour,Ref(T_period)),:]
# # required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_[in.(loads_.hour,Ref(T_period)),:], 0.05, 0)


G_N=100
ed_config = (
    remove_reserve_constraints = false,
)
G_MIP_GAP = 0.0001
G_RESERVE = 0.1

function test1(input_location = G_DEFAULT_LOCATION, reference_location = "./test/reference.csv")
    reference =  CSV.read(reference_location, DataFrame)
    
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_N, input_location)
    required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, G_RESERVE, 0)
    configs = generate_configurations(storage_df, required_reserve, required_energy_reserve, required_energy_reserve_cumulated)

    out = DataFrame(Dict(string(k) => solve_unit_commitment(
            gen_df,
            loads_multi_df,
            gen_variable_multi_df,
            G_MIP_GAP;
            v...).scalar[1, OBJECTIVE_VALUE] for (k,v) in configs
    ))
    # @infiltrate
    println(out == reference)
end

function test2()
    # configs_ = (base = configs[:base],)
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_N)
    required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, G_RESERVE, 0)
    configs = generate_configurations(storage_df, required_reserve, required_energy_reserve, required_energy_reserve_cumulated)

    out = [(solve_unit_commitment(
            gen_df,
            loads_multi_df,
            gen_variable_multi_df,
            G_MIP_GAP;
            v...).scalar[1, OBJECTIVE_VALUE],
            solve_economic_dispatch(
            gen_df,
            loads_multi_df,
            gen_variable_multi_df,
            G_MIP_GAP;
            merge(v, ed_config)...).scalar[1, OBJECTIVE_VALUE] ) for (k,v) in configs
    ]
    out = DataFrame(Config = collect(keys(configs)), UC= getindex.(out,1), EC = getindex.(out,2))
    out[!,:delta_percentual] .= (out.UC .- out.EC)./out.UC
    out[!,:delta_percentual_loq_mip_gap] .= out.delta_percentual .<=G_MIP_GAP
    println(out)
end

function test3()
    test1("./input/net_demand_base_case", "./test/reference_net_demand_base_case.csv")
end

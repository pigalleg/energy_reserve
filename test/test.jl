

include("../model/utils.jl")
include("../model/unit_commitment.jl")
include("../model/economic_dispatch.jl")
include("../notebooks/processing.jl")
# include("../plotting.jl")
using CSV
OBJECTIVE_VALUE = :objective_value
SCALAR = :scalar

G_N=100
ed_config = (
    remove_reserve_constraints = false,
)
G_MIP_GAP = 0.0001
G_RESERVE = 0.1

function test1(input_location = G_DEFAULT_LOCATION, reference_location = "./test/reference.csv")
    ref =  CSV.read(reference_location, DataFrame)
    
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_N, input_location)
    required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, G_RESERVE, 0)
    configs = generate_configurations(storage_df, required_reserve, required_energy_reserve, required_energy_reserve_cumulated)

    sol = DataFrame(Dict(string(k) => solve_unit_commitment(
            gen_df,
            loads_multi_df,
            gen_variable_multi_df,
            G_MIP_GAP;
            v...).scalar[1, OBJECTIVE_VALUE] for (k,v) in configs
    ))
    out = vcat(sol, ref, cols = :intersect)
    out[!,:header] = [:sol, :ref]
    out = permutedims(out, :header)
    out[!,:delta_percentual] .= (out.sol .- out.ref)./out.ref
    out[!,:delta_percentual_loq_mip_gap] .= out.delta_percentual .<=G_MIP_GAP
    println(all(out.delta_percentual_loq_mip_gap))
    return out
end

function test2()
    return test1("./input/net_demand_base_case", "./test/reference_net_demand_base_case.csv")
end

function test3(input_location = G_DEFAULT_LOCATION)
    # configs_ = (base = configs[:base],)
    gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(G_N, input_location)
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
    return out
end

function test4()
    return test3("./input/net_demand_case")
end

function test5(constrain_dispatch = false, reference_location = "./test/test5_reference.csv", variables_to_constrain = [:GEN])
    ref = change_type(change_type(CSV.read(reference_location, DataFrame), String, Symbol), String15, Symbol)
    days = [259]
    mip_gap = 0.000000001
    ed_config = Dict(:max_iterations => 10, :constrain_dispatch => constrain_dispatch, :remove_reserve_constraints => true, :variables_to_constrain => variables_to_constrain)

    μs = [(0,0), (0.1,0.1), (0.2,0.2), (0.3, 0.3), (0.4, 0.4), (0.5, 0.5), (0.6, 0.6), (0.7, 0.7), (0.75, 0.75), (0.8, 0.8), (0.85, 0.85), (0.9, 0.9), (0.95, 0.95), (1, 1)]
    configurations = [Symbol("base_ramp_storage_envelopes_up_$(replace(string(μ_up), "." => "_"))_dn_$(replace(string(μ_dn), "." => "_"))") for (μ_up, μ_dn) in μs]
    other_configs = [:base_ramp_storage_energy_reserve_cumulated]
    configurations = vcat(configurations, other_configs)

    s_ed = Dict()
    
    for day in days, k in configurations
        gen_df, loads_multi_df, gen_variable_multi_df, storage_df, random_loads_multi_df = generate_input_data(day, G_DEFAULT_LOCATION)
        required_reserve, required_energy_reserve, required_energy_reserve_cumulated = generate_reserves(loads_multi_df, gen_variable_multi_df, G_RESERVE)
        config = Dict(k => merge(v, ed_config) for (k,v) in generate_configuration(k, storage_df, required_reserve, required_energy_reserve, required_energy_reserve_cumulated))
        s_ed[(day,k)] = solve_economic_dispatch(
            gen_df,
            random_loads_multi_df,
            gen_variable_multi_df,
            mip_gap;
            config[k]...
        )
    end
    s_ed = merge_solutions(s_ed, [:day, :configuration])
    # CSV.write("test7_reference.csv", s_ed.scalar)
    out = leftjoin(
        s_ed.scalar[:,Not(:termination_status)], 
        rename(ref[:,Not(:termination_status)], :objective_value => :objective_value_ref),
        on = [:configuration, :day, :iteration])
    out[!,:delta_percentual] .= (out.objective_value .- out.objective_value_ref)./out.objective_value_ref
    out[!,:delta_percentual_loq_mip_gap] .= out.delta_percentual .<=mip_gap
    println(all(out.delta_percentual_loq_mip_gap))
    return out
end

function test6()
    return test5(true, "./test/test6_reference.csv", [:GEN])
end

function test7()
    return test5(true, "./test/test7_reference.csv", [:GEN, :CH, :DIS])
end
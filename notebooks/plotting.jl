using PlotlyJS

color_map=Dict(
    "battery" => "green",
    "solar_photovoltaic" => "gold",
    "net_generation" => "goldenrod",
    "onshore_wind_turbine" => "skyblue",
    "hydroelectric_pumped_storage" => "darkblue",
    "small_hydroelectric" => "cornflowerblue",
    "biomass" => "#6fc276",
    "natural_gas_fired_combined_cycle" => "grey",
    "natural_gas_fired_combustion_turbine" => "black",
    "total" => "purple",
    "required" => "#F0092",
    "envelope_down_MWh" => "blue",
    "envelope_up_MWh" => "red",
    "reserve_down_MW_eff" => "green",
    "reserve_up_MW_eff" => "green",
    "SOE_MWh" => "rgba(203, 213, 232, 1)",
)

color_discrete_map = (key) -> if haskey(color_map, key) color_map[key] else "red" end
map_stack_group = (x, exclude) -> 1*(x != exclude) + 2*(x == exclude)
TIMESTEP = :hour

function plot_fieldx_by_fieldy(df, fieldx, fieldy, title = nothing)
  out = plot(
    [scatter(
          x= df[df[!,fieldy] .== r, TIMESTEP], y=df[df[!,fieldy] .== r, fieldx],
          stackgroup="one", mode="lines", name = r,
          line=attr(width=1, color=color_discrete_map(r), shape = "line")
      ) for r in unique(df[!,fieldy])],
    Layout(yaxis_title="power MW", xaxis_title=string(TIMESTEP), title = string(fieldx))
  )
  if !isnothing(title)
      return relayout(out, title = title)
  end
  return out
end

function plot_reserve_by_fieldy(df, fieldx, fieldy, title = nothing)
  out = plot(
    [scatter(
          x= df[df[!,fieldy] .== r, TIMESTEP], y=df[df[!,fieldy] .== r, fieldx],
          stackgroup=map_stack_group(r,"required"), mode="lines", name = r,
          line=attr(width=1, color=color_discrete_map(r), shape = "line")
      ) for r in unique(df[!,fieldy])],
    Layout(yaxis_title="power MW", xaxis_title=string(TIMESTEP), title = string(fieldx))
  )
  if !isnothing(title)
      return relayout(out, title = title)
  end
  return out
end

function plot_supply_demand(supply, demand, title  = nothing)
  out = [plot_fieldx_by_fieldy(supply, :production_MW, :resource) plot_fieldx_by_fieldy(demand, :demand_MW, :resource)]
  if !isnothing(title)
      return relayout(out, title = title)
  end
  return out
end

function plot_reserve(reserve,  title  = nothing)
  out = [plot_reserve_by_fieldy(reserve, :reserve_up_MW, :resource) plot_reserve_by_fieldy(reserve, :reserve_down_MW, :resource)]
  if !isnothing(title)
      return relayout(out, title = title)
  end
  return out
end

function plot_battery_reserve_(battery_reserve, key)
  reserves = [scatter(
      x= battery_reserve.hour, y=battery_reserve[!,key],
      stackgroup=1, mode="lines", name = key,
      line=attr(width=1, shape = "line", color=color_discrete_map(string(key))), line_shape="vh"
  ) for key in [:SOE_MWh, key]]

  if :envelope_up_MWh in propertynames(battery_reserve)
      envelope = [scatter(
          x= battery_reserve.hour, y=battery_reserve[!,key],
          mode="lines", name = key,
          line=attr(width=1, shape = "line", color=color_discrete_map(string(key))), line_shape="vh"
          ) for key in [:envelope_up_MWh, :envelope_down_MWh]]
      union!(reserves, envelope)
  end
  
  return plot(reserves, Layout(yaxis_title="reserve MW", xaxis_title="hour", title = "All battery"))
end

function plot_battery_reserve(battery_reserve)
  return [plot_battery_reserve_(battery_reserve, :reserve_up_MW_eff) plot_battery_reserve_(battery_reserve, :reserve_down_MW_eff)]
end
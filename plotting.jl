using PlotlyJS

color_map=Dict(
    "battery" => "green",
    "solar_photovoltaic" => "gold",
    "onshore_wind_turbine" => "skyblue",
    "hydroelectric_pumped_storage" => "darkblue",
    "small_hydroelectric" => "cornflowerblue",
    "biomass" => "#6fc276",
    "natural_gas_fired_combined_cycle" => "grey",
    "natural_gas_fired_combustion_turbine" => "black",
    "total" => "purple",
    "required" => "#F0092",
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

function plot_battery_reserve(battery_reserve)
  reserve_keys = Dict(:up => :reserve_up_MW_eff, :down => :reserve_down_MW_eff)
  reserves_up = [scatter(
      x= battery_reserve.hour, y=battery_reserve[!,key],
      stackgroup=1, mode="lines", name = key,
      line=attr(width=1, shape = "line"), line_shape="vh"
  ) for key in [:SOE_MWh, reserve_keys[:up]]]

  reserves_down = [scatter(
      x= battery_reserve.hour, y=battery_reserve[!,key],
      stackgroup=1, mode="lines", name = key,
      line=attr(width=1, shape = "line"), line_shape="vh"
  ) for key in [:SOE_MWh, reserve_keys[:down]]]

  if :envelope_up_MWh in propertynames(battery_reserve)
      envelope = [scatter(
          x= battery_reserve.hour, y=battery_reserve[!,key],
          mode="lines", name = key,
          line=attr(width=1, shape = "line"), line_shape="vh"
          ) for key in [:envelope_up_MWh, :envelope_down_MWh]]
      union!(reserves_up, envelope)
      union!(reserves_down, envelope)
  end
  
  p1 = plot(reserves_up, Layout(yaxis_title="reserve MW", xaxis_title="hour", title = "All battery"))
  p2 = plot(reserves_down, Layout(yaxis_title="reserve MW", xaxis_title="hour", title = "All battery"))

  return [p1 p2]
end
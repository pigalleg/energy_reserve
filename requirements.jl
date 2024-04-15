using Pkg
dependencies = [
    "JuMP",
    "Gurobi",
    "DataFrames",
    "CSV",
    "Parquet2",
    "PlotlyJS",
    "Random",
    "Distributions",
    "MathOptInterface"
]

Pkg.add(dependencies)
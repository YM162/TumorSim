module TumorSim

#We export all the functions needed outside of the module
export create_scenario, create_treatment, simulate, launch_interactive_simulation, plot_genotypes, OncoSimulR_rfitness, get_TTP, get_diversity, Treatment, Scenario, launch_dashboard

#We import everything we need
using DrWatson
using Distributions
using StatsBase
using DataFrames
using Agents
using Random
using DataFrames
using DataStructures

#We include every file in the module.
include("Fitness/Fitness.jl")
include("Scenario/Scenario.jl")
include("Treatment/Treatment.jl")

include("TumorModel/TumorModel.jl")
include("Simulate/Simulate.jl")

include("Analysis/Analysis.jl")

include("Dashboard/Dashboard.jl")

end

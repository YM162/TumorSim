module TumorSim

#We export all the functions needed outside of the module
export create_scenario, create_treatment, simulate, launch_interactive_simulation, plot_genotypes, OncoSimulR_rfitness, get_TTP, get_diversity, Treatment, Scenario

#We import everything we need
using DrWatson
using Agents, Random
using Agents.DataFrames, Agents.Graphs
using Distributions
using InteractiveDynamics
using StatsBase
using ColorSchemes
using DataStructures
using DataFrames


#We include every file in the module.
include("Fitness/Fitness.jl")
#include("Fitness/OncoSimulR.jl")
include("Scenario/Scenario.jl")
include("Treatment/Treatment.jl")

include("TumorModel/TumorModel.jl")
include("Simulate/Simulate.jl")

include("Plotting/Plotting.jl")
include("Analysis/Analysis.jl")

end

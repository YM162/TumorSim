module TumorSim

#We export all the functions needed outside of the module
export create_scenario, create_treatment, simulate, get_TTP, get_diversity, get_resistant_fraction_on_detection, TreatmentObject, ScenarioObject, build_fitness_table

#We import everything we need
using DrWatson
using Agents, Random
using Agents.DataFrames, Agents.Graphs
using Distributions
using InteractiveDynamics
using StatsBase

#Mirar con el OhMyREPL el 
using DataStructures
using DataFrames

using TOML
function __init__()
    global Config = TOML.parse(open(projectdir("Config.toml")))
    nothing
end

#We include every file in the module.
include("Fitness/Fitness.jl")
include("Scenario/Scenario.jl")
include("Treatment/Treatment.jl")
include("TumorModel/TumorModel.jl")
include("Analysis/Analysis.jl")
include("Simulate/Simulate.jl")

end
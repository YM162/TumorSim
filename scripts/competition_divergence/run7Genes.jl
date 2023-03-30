using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

include(projectdir("scripts/competition_divergence/7Genes_NO_Restrictions.jl"))

include(projectdir("scripts/competition_divergence/7Genes_Example1.jl"))
include(projectdir("scripts/competition_divergence/7Genes_Example2.jl"))

include(projectdir("scripts/competition_divergence/7Genes_NK_Fitness_1.jl"))
include(projectdir("scripts/competition_divergence/7Genes_NK_Fitness_2.jl"))
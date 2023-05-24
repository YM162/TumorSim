TumorSim.jl
===========
![Build Status](https://github.com/YM162/TumorSim/actions/workflows/CI.yml/badge.svg)
[![codecov](https://codecov.io/gh/YM162/TumorSim/branch/main/graph/badge.svg?token=6YLRTP584L)](https://codecov.io/gh/YM162/TumorSim)

Julia code for the ***simulation of tumoral growth in multiple dimensions with a special focus on fitness landscapes and adaptive therapy.***

This code is used in Fontaneda, D. & Diaz-Uriarte, R. (2023). Adaptive therapy in cancer: the role of restrictions in the accumulation of mutations. BioRxiv, 2023.05.18(541330), doi: https://doi.org/10.1101/2023.05.18.541330


# How to install
<b>1.- Clone the git repo.</b>
```bash
git clone https://github.com/YM162/TumorSim.git
cd TumorSim
```
<b>2.- Install all required dependencies.</b>
First, you must install Julia: https://julialang.org/downloads/ . Then, install the dependencies:

```julia
using DrWatson
@quickactivate "TumorSim"

using Pkg; Pkg.instantiate()
```
<b>3.- (Optional) Build a sysimage for faster loading times.</b>
```bash
julia scripts/build.jl
```
# How to use
<b>1.- (Optional) Launch the sysimage and set the number of cores for paralellization.</b>
```bash
julia -q --sysimage build/TumorSim.so -p 4
```
<b>2.- Import the DrWatson and TumorSim packages.</b>
```julia
using DrWatson
@quickactivate "TumorSim"
using TumorSim
```
<b>3.- Create a scenario defining the starting state of the simulation.</b>
```julia
# create_scenario(dims, Initial Cells, cell_positions).
# Cell position can be "center", "random" or an array specifying the exact positions.

scenario_1D = create_scenario((1000000,),10,"center")

scenario_2D = create_scenario((1000,1000),10,"center")

scenario_3D = create_scenario((100,100,100),10,"center")
```
<b>4.- Create a fitness landscape.</b> 
```julia
#Must be a dictionary with genotype keys mapping to mitosis probabilities. If a genotype is not in the dictionary, it is asumed to have a mitosis probability of 0.

fitness = Dict([0,0,0]=>0.027,
                [1,0,0]=>0.030,
                [0,1,0]=>0.033,
                [1,1,0]=>0.046,
                [1,1,1]=>0.025)

#We can also use the OncoSimulR rfitness function. 
#See ./scripts/AT_in_cancer_the_role_of_restrictions_2023/Simulate_7Genes_NK_Fitness_1.jl for a working example.
```
>For more information about generating random fitness landscapes take a look at the [OncoSimulR documentation](https://www.bioconductor.org/packages/release/bioc/vignettes/OncoSimulR/inst/doc/OncoSimulR.html#9_Generating_random_fitness_landscapes).

<b>5.- Create a Treatment.</b>
```julia
#Treatment(Detection size, Starting size (% of detection size), Pausing size (% of detection size), Gene of resistance, kill_rate)

adaptive_therapy = create_treatment(10000, 1.0, 0.5, 3, 0.75) 

continuous_therapy = create_treatment(10000, 1.0, 0, 3, 0.75) 

#After the Detection size has been reached, treatment cycles starting and pausing at the specified sizes will begin, killing kill_rate% of the susceptible cells that try to reproduce.
```

<b>8.- Create a vector of all the scenarios you want to test.</b>
```julia
#We need to specify death rate (% of WT mitosis probability), mutation rate, starting scenario, fitness landscape, treatment, migration rate and interaction rule (:contect or :hierarchical_voter)

parameters = Dict(
    "death_rate" => [0,0.3],
    "mutation_rate" => [0.01,0.015,0.02],   
    "scenario" => [scenario_2D,scenario_3D], 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy,continuous_therapy],
    "migration_rate" => [0.1],
    "interaction_rule" => [:contact],
    "seed" => map(abs,rand(Int64,1000))
)

#DrWatson's dict_list allows us to expand the vectors in our original dict to produce a list with all of the possible combinations of parameters.

parameter_combinations = dict_list(parameters)
```
<b>9a.- Run a single simulation.</b>
```julia
#We specify the number of steps each simulation will go through. The simulation will stop early if all cells die or if we reach 1.5 * treatment.detection_size (resistance was aquired).
steps=5000

#Get a dictionary with the results of the simulation.
result = simulate(parameter_combinations[1],steps)

#We save it to disk in the BSON format
using BSON
filepath = datadir("simulations","example_sim.bson")
bson(filepath,Dict("result" => result))

```
<b>9b.- Run all the simulations at once using parallelization.</b>
```julia
using Distributed
using ProgressMeter

#We use progress_pmap from the ProgressMeter package to take advantage of parallelization and perform the simulations.
p = Progress(length(parameter_combinations), 
            barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)

results = progress_pmap(simulate,
                        parameter_combinations,
                        fill(steps,length(parameter_combinations)),progress=p)


#We again save all of the simulations in a dataframe in the BSON format.
using BSON
using DataFrames

df = DataFrame(results)
filepath = datadir("simulations","example_multi_sim.bson")
bson(filepath,Dict("df" => df))
```


# License
All code is distributed under the GNU GPL v. 3.0 license.

# Citation
A pre-print of our manuscript using this model can be found in bioRxiv:
* Fontaneda, D. & Diaz-Uriarte, R. (2023). Adaptive therapy in cancer: the role of restrictions in the accumulation of mutations. BioRxiv, 2023.05.18(541330), doi: https://doi.org/10.1101/2023.05.18.541330


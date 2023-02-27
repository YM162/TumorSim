TumorSim.jl
===========
![Build Status](https://github.com/YM162/TumorSim/actions/workflows/CI.yml/badge.svg)
[![codecov](https://codecov.io/gh/YM162/TumorSim/branch/main/graph/badge.svg?token=6YLRTP584L)](https://codecov.io/gh/YM162/TumorSim)

***Simulation of tumoral growth in multiple dimensions with a special focus on fitness landscapes and adaptive therapy.***

# How to install
<b>1.- Clone the git repo.</b>
```bash
git clone https://github.com/YM162/TumorSim.git
cd TumorSim
```
<b>2.- (Optional) Build a sysimage for faster loading times.</b>
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

scenario_0D = create_scenario(1000000,10)

scenario_1D = create_scenario((1000000,),10,"center")

scenario_2D = create_scenario((1000,1000),10,"center")

scenario_3D = create_scenario((100,100,100),10,"center")
```
<b>4.- Create a fitness landscape.</b> 
```julia
#Must be a dictionary with genotype keys mapping to fitness values.

fitness=Dict([0,0,0]=>1, 
             [1,0,0]=>1.3,
             [0,1,0]=>1.2,
             [1,1,0]=>1.5,
             [1,1,1]=>1.2)

#We can also use the OncoSimulR rfitness function.
fitness=OncoSimulR_rfitness(g=3,c=1,sd=0.5)
```
>For more information about generating random fitness landscapes take a look at the [OncoSimulR documentation](https://www.bioconductor.org/packages/release/bioc/vignettes/OncoSimulR/inst/doc/OncoSimulR.html#9_Generating_random_fitness_landscapes).

<b>5.- Create a Treatment.</b>
```julia
#Treatment(Detection size, Starting size, Pausing size, Gene of resistance, kill_rate)

adaptive_therapy = create_treatment(3000, 2000, 1000, 3, 0.75) 

continuous_therapy = create_treatment(3000, 2000, 0, 3, 0.75) 

#After the Detection size has been reached, treatment cycles starting and pausing at the specified sizes will begin, killing kill_rate% of the susceptible cells that try to reproduce.
```

<b>8.- Create a vector of all the scenarios you want to test.</b>
```julia
#We need to specify pr (Base proliferation rate), dr (Death rate) and mr (Mutation rate)

parameters = Dict(
    "pr" => [0.01,0.02,0.03,0,04,0.05],
    "dr" => [0,0.2,0.4,0.6,0.8],
    "mr" => [0.001,0.005,0.01,0.05,0.1],   
    "scenario" => [scenario_0D,scenario_1D,scenario_2D,scenario_3D], 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy,continuous_therapy],
    "seed" => map(abs,rand(Int64,1000))
)

#DrWatson's dict_list allows us to expand the vectors in our original dict to produce a list with all of the possible combinations of parameters.

parameter_combinations = dict_list(parameters)
```
<b>9a.- Run a single simulation.</b>
```julia
#We specify the number of steps each simulation will go through. The simulation will stop early if all cells die or if we reach 1.5 * treatment.detection_size (resistance was aquired).
steps=1000

#Get a dictionary with the results of the simulation.
result = simulate(parameter_combinations[1],steps)

#We save it to disk using DrWatson's safesave and savename functions.
safesave(datadir("simulations", savename(result, "jld2")),result)

```
<b>9b.- Run all the simulations at once using parallelization.</b>
```julia
using Distributed
using ProgressMeter
#We use pmap to take advantage of parallelization to perform the simulations.
results = @showprogress pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)))

#We again use DrWatson's safesave and savename functions to save each simulation to disk.
for (i, d) in enumerate(parameter_combinations)
    safesave(datadir("simulations", savename(results[i], "jld2")), results[i])
end
```
<b>10.- Collect the data from all the simulations you saved.</b>
```julia
df = collect_results(datadir("simulations"))
```
# <b>tumor-sim</b> </br>
Simulation of tumoral growth using Agents.jl </br>

# Current features:
* Scenarios for 0D, 1D, 2D and 3D.
* Continuous and adaptive therapy.
* Fitness landscapes with arbitrary number of genes.

# How to install
```julia
using Pkg;Pkg.add(url="https://github.com/YM162/TumorSim.git")
```
# How to use:
(Optional) Build a sysimage and launch a julia prompt with it. Enable paralellyzation with -p n_cores
```
julia build.jl
julia -q -p 4 --sysimage build/TumorSim.so
```
1.- Import the TumorSim package
```julia
using TumorSim
```
2.- Create a scenario with the dimensions, initial number of cells and initial position of cells.
```julia
scenario = create_scenario((100,100,100),10,"center")
```
3.- Create a fitness landscape
```julia
fitness=Dict([0,0,0]=>1, 
             [1,0,0]=>1.3,
             [0,1,0]=>1.2,
             [1,1,0]=>1.5,
             [1,1,1]=>1.2)
```
4.- Create a Treatment specifying detection size, starting size, pausing size, gene of resistance and kill rate.
```julia
adaptive_therapy = Treatment(3000, 2000, 1000, 3, 0.75, false, false) 
```
5.- Define the data to collect, in this case, the cell count for all the genotypes in fitness.
```julia
agent_collect = [(:genotype, f) for f in genotype_fraction_function_generator(fitness)]
```
6.- Initiallize the model
```julia
model = model_init(pr=0.027, #Proliferation rate
                    dr=0.015, #Death rate
                    mr=0.01, #Mutation rate
                    scenario=scenario, 
                    fitness=fitness, 
                    treatment=adaptive_therapy, 
                    seed=0)

```
7.- Run the model
```julia
adata, _ = run!(model, agent_step!, model_step!, steps; adata = agent_collect)
```
8.- Clean up the dataframe
```julia
genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]

pushfirst!(genotypes,"step")
rename!(adata,genotypes)
```
#Copy from the code of tumor-sim.ipynb as it was on 17/11/22

#We import the common functions and model.
include("./functions.jl")
using GLMakie
GLMakie.activate!() 
#Number of genes of each cell. We need to define it earlier in order to collect the data. Can go up to 10000 easily AS LONG as we dont measure in every timestep.
ngenes=3 

#
# Ejemplo AND con 3 genes:
#
#  WT             000:WT  ->1
#  /\             001:C   ->0
# A  B            010:B   ->1.2
#  \/ (and)       011:BC  ->0
#   C             100:A   ->1.3
#                 101:AC  ->0
#                 110:AB  ->1.5
#                 111:ABC ->1.2
#

fitness=Dict(0=>1, #Fitness of each of the genotypes, the key is the binary representation of the genotype in decimal form.
            1=>0, #Fitness is expressed as a multiplicative effect on pr, but this can be changed.
            2=>1.2,
            3=>0,
            4=>1.3,
            5=>0,
            6=>1.5,
            7=>1.2)


#OncoSimulR integration! We can use the rfitness function to generate random fitness landscapes
#fitness=OncoSimulR_rfitness(g=ngenes,c=0.5,sd=1) 

#If a genotype fitness is not specified, it can have a default value (0 for not viable or 1 for same effect as WT)
fitness=DefaultDict(1,fitness) 

#We load a scenario from an image
h,w,cell_pos,wall_pos = get_scenario("./scenarios/example2.bmp") #I have added some examples (1 to 4) and some templates.

#We define a Treatment
#I need to include the posibility of making it time-dependent instead of tumor size-dependent.
adaptive_therapy = Treatment(1000, #Pausing size: treatment will stop when this cell_count is reached
                            2000, #Starting size: treatment will start when this cell_count is reached
                            3, #Gene of resistance
                            0.75, #Kill rate
                            false) #Initial state of the treatment

                            
continuous_therapy = Treatment(0, #Pausing size
                            2000, #Starting size
                            3, #Gene of resistance
                            0.75, #Kill rate
                            false) #Initial state of the treatment

#We initialize the model.
model = model_init(pr=0.027, #Proliferation rate
                    dr=0.015, #Death rate
                    mr=0.005, #Mutation rate
                    h=h, #Height of the grid
                    w=w, #Witdh of the grid
                    cell_pos=cell_pos, #The initial positions of the cells
                    wall_pos=wall_pos, #Places the cells can not go
                    ngenes=ngenes, #The number of genes of each cell
                    fitness=fitness, #The fitness of each genotype
                    treatment=adaptive_therapy, #The treatment we are going to use for the simulation
                    seed=3) #Seed to get reproducible results

println(model)

#we collect the number of cells of each genotype that are alive
to_collect = [(:genotype, f) for f in genotype_fraction_function_generator(ngenes)]

#we run the simulation
steps= 1700

data, _ = run!(model, agent_step!, model_step!, steps; adata = to_collect)

#we rename the columns to get a clean "data" DataFrame.
genotypes = [string(reverse(digits(i, base=2, pad=ngenes))) for i in 0:((2^ngenes)-1)]
pushfirst!(genotypes,"step")
rename!(data,genotypes)


#We plot the genotype of each cell with a different color.

#Right now we plot the walls as a heatarray because its te easiest thing i found to just get it working.
heatarray = :wall_matrix
heatkwargs = (colorrange = (0, 1), colormap = :grayC)

#Functions to get a different color for each genotype
genotypecolor(a) = get(colorschemes[:hsv], bit_2_int(a.genotype), (0,(2^ngenes)+1))
genotypecolor_legend(a) = get(colorschemes[:hsv], a, (0,(2^ngenes)+1))

#We make a static plot
figure, _ = abmplot(model;ac = genotypecolor,as=8,am='■',heatarray,heatkwargs)

#We create a legend for the genotypes
genotypes = [filter(x -> !isspace(x), string(reverse(digits(i, base=2, pad=ngenes)))) for i in 0:((2^ngenes)-1)]
Legend(figure[1, 2],
    [MarkerElement(color = genotypecolor_legend(a), marker = '■', markersize = 15, strokecolor = :black) for a in 0:(2^ngenes)-1],
    genotypes,
    patchsize = (20, 20), rowgap = 1)


#We display the figure
display(figure)

#And lastly we can make plots of both the total number of cells of each genotype
genotypes = [string(reverse(digits(i, base=2, pad=ngenes))) for i in 0:((2^ngenes)-1)]
stacked = stack(data,genotypes)

stacked |>
@vlplot(:area, x=:step, y={:value, stack=:zero}, color="variable:n")

#And the relative number of cells of each genotype
#stacked |>
#@vlplot(:area, x=:step, y={:value, stack=:normalize}, color="variable:n")
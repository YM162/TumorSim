#We import the common functions and model.
println("""     ---- INTERACTIVE DEMO ----

* Space: 2D Grid (100x100)
* Initial scenario: scenarios/example3.bmp
* Fitness landscape (0 if not specified): 
    WT -> 1
    A -> 1.3
    B -> 1.2
    AB -> 1.5
    ABC -> 1.2
* Treatment: Adaptive therapy (start: 2000 cells, pause: 1000 cells, kill_rate: 0.75, gene of resistance: C)
""")
println("LOADING LIBRARIES: Warmup may take up to 2 minutes because Julia. Program will launch eventually in a new window")
include("./functions.jl")
#We use GLMakie to use the interactive plot
using GLMakie
GLMakie.activate!() 

function launch_from_scenario(scenario_path)
    #We load a scenario from an image
    h,w,cell_pos,wall_pos = get_scenario(scenario_path) #I have added some examples (1 to 4) and some templates.

    size_marker = 800/max(h,w)
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
    #                 111:ABC ->2
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
    fitness=DefaultDict(0,fitness) 

   
    #We define a Treatment
    adaptive_therapy = Treatment(1000, #Pausing size
                                2000, #Starting size
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
    println("Model:")
    println(model)

    #we collect the number of cells of each genotype that are alive
    to_collect = [(:genotype, f) for f in genotype_fraction_function_generator(ngenes)]

    #we run the simulation
    steps= 0

    data, _ = run!(model, agent_step!, model_step!, steps; adata = to_collect)

    #we rename the columns to get a clean "data" DataFrame.
    genotypes = [string(reverse(digits(i, base=2, pad=ngenes))) for i in 0:((2^ngenes)-1)]
    pushfirst!(genotypes,"step")
    rename!(data,genotypes)

    #We plot the genotype of each cell with a different color.

    #Functions to get a different color for each genotype
    genotypecolor(a) = get(colorschemes[:hsv], bit_2_int(a.genotype), (0,(2^ngenes)+1))
    genotypecolor_legend(a) = get(colorschemes[:hsv], a, (0,(2^ngenes)+1))

    #Right now we plot the walls as a heatarray because its te easiest thing i found to just get it working.
    heatarray = :wall_matrix
    heatkwargs = (colorrange = (0, 1), colormap = :grayC)

    #We create the sliders to modify the model
    params = Dict(
        :pr => 0:0.001:0.1,
        :dr => 0:0.001:0.1,
        :mr => 0:0.001:0.05
    )

    #We make an interactive plot
    figure, ax, abmobs = abmplot(model;agent_step! = agent_step!, model_step! = model_step!, ac = genotypecolor,as=size_marker,am='■',heatarray,heatkwargs,params)

    #We create a legend for the genotypes
    genotypes = [filter(x -> !isspace(x), string(reverse(digits(i, base=2, pad=ngenes)))) for i in 0:((2^ngenes)-1)]
    Legend(figure[1, 2],
        [MarkerElement(color = genotypecolor_legend(a), marker = '■', markersize = 15, strokecolor = :black) for a in 0:(2^ngenes)-1],
        genotypes,
        patchsize = (20, 20), rowgap = 1)

    #We display the figure in a new window
    display(figure)
end

function run_interactive()
    #scenario = ones(RGB{N0f8},h,w)
    #We load scenario 3
    scenario = FileIO.load("./scenarios/example2.bmp")
    scenario = rotr90(scenario)

    reset_scenario = copy(scenario)
    x = Observable(scenario)
    fig2, ax2, plotobj2 = image(x,interpolate=false)
    ##
    pcolor = Observable(RGB{N0f8}(1,1,1))
    ax2.title="Change the starting scenario using the tools below."
    fig2[2, 1] = buttongrid = GridLayout(tellwidth = false)

    buttons = buttongrid[1, 1:5] = [Button(fig2, label = "Reset figure"),Button(fig2, label = "Draw Walls"),Button(fig2, label = "Draw cells"),Button(fig2, label = "Eraser"),Button(fig2, label = "Start simulation")]

    ##
    display(fig2)
    ##
    on(buttons[1].clicks) do n
        scenario = copy(reset_scenario)
        x[] = scenario
    end
    on(buttons[2].clicks) do n
        pcolor[] = RGB{N0f8}(0,0,0)
        ax2.title="Tool selected: Wall creator"
    end
    on(buttons[3].clicks) do n
        pcolor[] = RGB{N0f8}(0.133,0.694,0.298)
        ax2.title="Tool selected: Cell creator"
    end
    on(buttons[4].clicks) do n
        pcolor[] = RGB{N0f8}(1,1,1)
        ax2.title="Tool selected: Eraser"
    end
    on(buttons[5].clicks) do n #this one need to save the image.
        ax2.title="LAUNCHING SIMULATION: Please Wait"
        FileIO.save("tmp/last_scenario.bmp",rotl90(x[]))
        launch_from_scenario("tmp/last_scenario.bmp")
    end
    ##
    on(events(fig2).mousebutton, priority=0) do event
        if event.button == Mouse.left
            scenario = x[]
            xpos, ypos = mouseposition(ax2.scene)
            scenario[floor(Int,xpos+1),floor(Int,ypos+1)] = pcolor[]
            x[] = scenario
        end
    end
end

run_interactive()

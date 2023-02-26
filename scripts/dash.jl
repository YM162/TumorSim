using DrWatson
@quickactivate "TumorSim"
using Blink
using InteractiveDynamics
using Agents
using TumorSim

w = Window(Dict("minWidth"=> 800,"minHeight"=> 450))

@js_ w document.open()
@js_ w document.write($(String(read(srcdir("Dashboard/dashhtml/index.html")))))
load!(w,srcdir("Dashboard/dashhtml/assets/css/style.css"))
@js_ w document.write($(String(read(srcdir("Dashboard/dashhtml/singlethree.html")))))


function get_abmobs(;fitness=Dict([0,0,0]=>0.027, [1,0,0]=>0.031,[0,1,0]=>0.035,[1,1,0]=>0.040,[1,1,1]=>0.040),
    scenario = create_scenario((100,100,100),10,"center",false),
    treatment = create_treatment(3000, 1, 0.5, 3, 0.75),
    cost_of_resistance=0.1,
    death_rate=0.1,
    mutation_rate=0.01,
    migration_rate=0.01,
    seed=0) 

    model = TumorSim.TumorModel.model_init(cost_of_resistance=cost_of_resistance, death_rate=death_rate, mutation_rate=mutation_rate, migration_rate=migration_rate,
        scenario=scenario, fitness=fitness,treatment=treatment, 
        seed=seed)

    abmobs = ABMObservable(model;agent_step! = TumorSim.TumorModel.agent_step!,model_step! = TumorSim.TumorModel.model_step!)
    return abmobs
end
abmobs = get_abmobs()
resetargs = Any[0.2, 0.01, 0.01, 485542213, Any[100, 100, 100], 100, "center", "false", 3000, 1, 0.5, 0.75, "4", 0.2]
handle(w, "changepage") do args
    
    global singlesetup
    @js_ w document.open()
    
    @js_ w document.write($(String(read(srcdir("Dashboard/dashhtml/"*args*".html")))))
    
    sleep(0.5)
    load!(w,srcdir("Dashboard/dashhtml/assets/css/style.css"))
    if args == "singlesimulation"
        
        println("Setting up single simulation")
        @js_ w document.write($(String(read(srcdir("Dashboard/dashhtml/singlethree.html"))))) #Double sending magic makes it work.

        @js_ w document.write($(String(read(srcdir("Dashboard/dashhtml/singlethree.html")))))
        sleep(0.5)
        @js_ w generate_plot1()
        
        @js_ w window.initvox()
        sleep(0.5)
        global abmobs = get_abmobs()
        step_and_plot(abmobs,0)
        
    end
    
end

#For single simulations

function step_and_plot(abmobs,s)
    step!(abmobs,s)
    updategraph(abmobs.model[])
end

oldpositions = [[[1,1,1],0]] #Placeholder for first iteration
function updategraph(model)
    global oldpositions
    genotype_bits = [x for x in sort!([x for x in keys(model.fitness)],by=x -> TumorSim.Fitness.bit_2_int(BitArray(x)))]
    order = Dict(zip(genotype_bits,1:length(genotype_bits)))
    
    positions = [[[a for a in x.pos],order[x.genotype]] for x in collect(allagents(model))]

    changes=vcat(setdiff(positions,oldpositions),[[x,0] for x in [x[1] for x in oldpositions] if x ∉ [x[1] for x in positions]])

    @js_ w window.world.setVoxelArray($(changes));
    global oldpositions = positions
    @js_ w window.updateVoxelGeometry(1,1,1)
    @js_ w window.requestRenderIfNotRequested()

    return nothing
end

function create_single(args)
    death_rate = Float64(args[1])
    mutation_rate = Float64(args[2])
    migration_rate = Float64(args[3])
    seed = args[4]
    size = args[5]
    initial_cells = args[6]
    initial_positions = args[7]
    mix = args[8]
    detecting_size = args[9]
    starting_size = Float64(args[10])
    pausing_size = Float64(args[11])
    kill_rate = Float64(args[12])
    restriction_scenario = args[13]
    cost_of_resistance = Float64(args[14])

    if restriction_scenario == "1"
        fitness = Dict([0,0,0]=>0.027, [1,0,0]=>0.031,[0,1,0]=>0.035,[0,0,1]=>0.027,[1,1,0]=>0.040,[0,1,1]=>0.035,[1,0,1]=>0.031,[1,1,1]=>0.040)
    end
    if restriction_scenario == "2"
        fitness = Dict([0,0,0]=>0.027, [1,0,0]=>0.031,[0,1,0]=>0.035,[1,1,0]=>0.040,[0,1,1]=>0.035,[1,0,1]=>0.031,[1,1,1]=>0.040)
    end
    if restriction_scenario == "3"
        fitness = Dict([0,0,0]=>0.027, [1,0,0]=>0.031,[0,1,0]=>0.035,[1,1,0]=>0.040,[0,1,1]=>0.035,[1,1,1]=>0.040)
    end
    if restriction_scenario == "4"
        fitness = Dict([0,0,0]=>0.027, [1,0,0]=>0.031,[0,1,0]=>0.035,[1,1,0]=>0.040,[1,1,1]=>0.040)
    end
    if restriction_scenario == "5"
        fitness = Dict([0,0,0]=>0.027, [1,0,0]=>0.031,[1,1,0]=>0.040,[1,1,1]=>0.040)
    end

    abmobs = get_abmobs(;fitness=fitness,
                        scenario = create_scenario(Tuple(size),initial_cells,initial_positions,parse(Bool,mix)),
                        treatment = create_treatment(detecting_size, starting_size, pausing_size, 3, kill_rate),
                        cost_of_resistance=cost_of_resistance,
                        death_rate=death_rate,
                        mutation_rate=mutation_rate,
                        migration_rate=migration_rate,
                        seed=seed)
    return abmobs
end

handle(w, "create_single") do args
    global abmobs
    @show args
    abmobs = create_single(args)
    
    global resetargs = copy(args)
    step_and_plot(abmobs,0)
end


handle(w, "resetabmobs") do args
    global resetargs
    global abmobs = create_single(resetargs)
    step_and_plot(abmobs,0)
end

handle(w, "nrun") do args
    step_and_plot(abmobs,1)
end

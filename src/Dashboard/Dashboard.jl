module Dashboard
    export launch_dashboard

    using DrWatson
    using Blink
    using InteractiveDynamics
    using Agents
    using TumorSim
    using DataFrames
    mutable struct DashboardStruct
        w::Window
        single_abmobs::Any
        single_oldpositions::Vector
        single_resetargs::Vector
    end

    function get_abmobs(;fitness=Dict([0,0,0]=>0.027, [1,0,0]=>0.031,[0,1,0]=>0.035,[1,1,0]=>0.040,[1,1,1]=>0.040),
        scenario = create_scenario((100,100,100),100,"center",false),
        treatment = create_treatment(3000, 1, 0.5, 3, 0.75),
        cost_of_resistance=0.2,
        death_rate=0.1,
        mutation_rate=0.01,
        migration_rate=0.01,
        seed=0) 
    
        model = TumorSim.TumorModel.model_init(cost_of_resistance=cost_of_resistance, death_rate=death_rate, mutation_rate=mutation_rate, migration_rate=migration_rate,
            scenario=scenario, fitness=fitness,treatment=treatment, 
            seed=seed)
        
        agent_collect::Array = [(:genotype, f) for f in TumorSim.genotype_fraction_function_generator(fitness)]
        abmobs = ABMObservable(model;agent_step! = TumorSim.TumorModel.agent_step!,model_step! = TumorSim.TumorModel.model_step!, adata = agent_collect)
        return abmobs
    end

    function step_and_plot(abmobs,s)
        step!(abmobs,s)
        updategraph(abmobs.model[])
    end

    function updategraph(model)
        genotype_bits = [x for x in sort!([x for x in keys(model.fitness)],by=x -> TumorSim.Fitness.bit_2_int(BitArray(x)))]
        order = Dict(zip(genotype_bits,1:length(genotype_bits)))
        
        positions = [[[a for a in x.pos],TumorSim.bit_2_int(x.genotype)+1] for x in collect(allagents(model))]
    
        changes=vcat(setdiff(positions,TumorSim_Dashboard.single_oldpositions),[[x,0] for x in [x[1] for x in TumorSim_Dashboard.single_oldpositions] if x ∉ [x[1] for x in positions]])
    
        @js_ TumorSim_Dashboard.w window.world.setVoxelArray($(changes));
        TumorSim_Dashboard.single_oldpositions = positions
        @js_ TumorSim_Dashboard.w window.updateVoxelGeometry(1,1,1)
        @js_ TumorSim_Dashboard.w window.requestRenderIfNotRequested()
    
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
        #Ahora está hard-coded, hay que pensar una manera de permitir introducir otros escenarios al usuario
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
    
    function clearvox()
        changes = [[x,0] for x in [x[1] for x in TumorSim_Dashboard.single_oldpositions]]
        @js_ TumorSim_Dashboard.w window.world.setVoxelArray($(changes));
        @js_ TumorSim_Dashboard.w window.updateVoxelGeometry(1,1,1)
        @js_ TumorSim_Dashboard.w window.requestRenderIfNotRequested()
        TumorSim_Dashboard.single_oldpositions = [[[1,1,1],0]]
    end

    function launch_dashboard(test=false)
        w = Window(Dict("minWidth"=> 800,"minHeight"=> 450,"width"=> 1250,"height"=> 675))

        global TumorSim_Dashboard = DashboardStruct(w,
                                            get_abmobs(),
                                            [[[1,1,1],0]],
                                            [0.2, 0.01, 0.01, 485542213, Any[100, 100, 100], 100, "center", "false", 3000, 1, 0.5, 0.75, "4", 0.2])


        @js_ TumorSim_Dashboard.w document.open()
        @js_ TumorSim_Dashboard.w document.write($(String(read(srcdir("Dashboard/dashhtml/index.html")))))
        load!(TumorSim_Dashboard.w,srcdir("Dashboard/dashhtml/assets/css/style.css"))
        @js_ TumorSim_Dashboard.w document.write($(String(read(srcdir("Dashboard/dashhtml/singlethree.html")))))
        

        handle(TumorSim_Dashboard.w, "changepage") do args
            
            @js_ TumorSim_Dashboard.w document.open()
            @js_ TumorSim_Dashboard.w document.write($(String(read(srcdir("Dashboard/dashhtml/"*args*".html")))))
            load!(TumorSim_Dashboard.w,srcdir("Dashboard/dashhtml/assets/css/style.css"))

            if args == "singlesimulation"
                TumorSim_Dashboard.single_abmobs = get_abmobs()
                clearvox()
                println("Setting up single simulation")
                @js_ TumorSim_Dashboard.w document.write($(String(read(srcdir("Dashboard/dashhtml/singlethree.html"))))) #Double sending magic makes it work.
        
                @js_ TumorSim_Dashboard.w document.write($(String(read(srcdir("Dashboard/dashhtml/singlethree.html")))))
                sleep(0.5)
                js(TumorSim_Dashboard.w,Blink.JSString("""generate_plot1($(TumorSim_Dashboard.single_abmobs.model[].ngenes),$(TumorSim_Dashboard.single_abmobs.model[].treatment.detecting_size),$(Int(floor(TumorSim_Dashboard.single_abmobs.model[].treatment.detecting_size*TumorSim_Dashboard.single_abmobs.model[].treatment.starting_size))),$(Int(floor(TumorSim_Dashboard.single_abmobs.model[].treatment.detecting_size*TumorSim_Dashboard.single_abmobs.model[].treatment.starting_size*TumorSim_Dashboard.single_abmobs.model[].treatment.pausing_size))))"""),callback=false)
                
                @js_ TumorSim_Dashboard.w window.initvox()
                sleep(0.5)
                
                step_and_plot(TumorSim_Dashboard.single_abmobs,0)
                
            end
            
        end
   
        handle(w, "create_single") do args
            @show args
            TumorSim_Dashboard.single_abmobs = create_single(args)
            clearvox()
            TumorSim_Dashboard.single_resetargs = copy(args)
            step_and_plot(TumorSim_Dashboard.single_abmobs,0)
            js(TumorSim_Dashboard.w,Blink.JSString("""generate_plot1($(TumorSim_Dashboard.single_abmobs.model[].ngenes),$(TumorSim_Dashboard.single_abmobs.model[].treatment.detecting_size),$(Int(floor(TumorSim_Dashboard.single_abmobs.model[].treatment.detecting_size*TumorSim_Dashboard.single_abmobs.model[].treatment.starting_size))),$(Int(floor(TumorSim_Dashboard.single_abmobs.model[].treatment.detecting_size*TumorSim_Dashboard.single_abmobs.model[].treatment.starting_size*TumorSim_Dashboard.single_abmobs.model[].treatment.pausing_size))))"""),callback=false)
        end
        
        handle(w, "resetabmobs") do args
            TumorSim_Dashboard.single_abmobs = create_single(TumorSim_Dashboard.single_resetargs)
            clearvox()
            step_and_plot(TumorSim_Dashboard.single_abmobs,0)
            js(TumorSim_Dashboard.w,Blink.JSString("""generate_plot1($(TumorSim_Dashboard.single_abmobs.model[].ngenes),$(TumorSim_Dashboard.single_abmobs.model[].treatment.detecting_size),$(Int(floor(TumorSim_Dashboard.single_abmobs.model[].treatment.detecting_size*TumorSim_Dashboard.single_abmobs.model[].treatment.starting_size))),$(Int(floor(TumorSim_Dashboard.single_abmobs.model[].treatment.detecting_size*TumorSim_Dashboard.single_abmobs.model[].treatment.starting_size*TumorSim_Dashboard.single_abmobs.model[].treatment.pausing_size))))"""),callback=false)
        end
        
        handle(w, "nrun") do args
            step_and_plot(TumorSim_Dashboard.single_abmobs,1)
            
            adata = TumorSim_Dashboard.single_abmobs.adf[]
            genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(TumorSim_Dashboard.single_abmobs.model[].fitness)],by=x -> TumorSim.bit_2_int(BitArray(x)))]
            lastdata = [adata[end,i] for i in 1:ncol(adata)]
            for gen in 1:(length(lastdata)-1)
                id = TumorSim.bit_2_int(BitArray(eval(Meta.parse(genotypes[gen]))))
                number = lastdata[gen+1]
                js(TumorSim_Dashboard.w,Blink.JSString("""Plotly.extendTraces("plot1",{y:[[$number]]},[$id])"""),callback=false)
            end
            js(TumorSim_Dashboard.w,Blink.JSString("""window.inprocess=false"""),callback=false)
        end
        
    end
end
using .Dashboard
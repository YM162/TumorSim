


@reactive mutable struct LaunchSimulationsPage <: ReactiveModel
    
    tableData::R{DataTable} = DataTable(DataFrame(Name=[],Progress=[]))
    credit_data_pagination::DataTablePagination = DataTablePagination(rows_per_page=10)

    name::R{String} = ""

    simulate_button::R{Int} = 0
    cancel_button::R{Int} = 0

    predicted_simulations::R{Int} = 20

    scenario_list::R{Vector{Int64}} = [0, 1, 2, 3]
    scenario::R{Vector{Int64}} = [3]
                                          
    detecting_size::R{RangeData{Int64}} = RangeData(3000:3000)
    starting_size::R{RangeData{Int64}} = RangeData(800:800)
    pausing_size::R{RangeData{Int64}} = RangeData(500:500)

    death_rate::R{RangeData{Int64}} = RangeData(500:500)
    mutation_rate::R{RangeData{Int64}} = RangeData(10:10)

    cost_of_resistance::R{RangeData{Int64}} = RangeData(200:200)

    kill_rate::R{RangeData{Int64}} = RangeData(750:750)

    repetitions::R{Int64} = 10

end

function launch_simulations_ui(model::LaunchSimulationsPage)
    
    @async while true
        
        simulations = readdir(projectdir("logs","progress"))
        df = DataFrame(Name=[],Progress=[])
        for simulation in simulations
            append!(df, DataFrame(Name=simulation,Progress=split(readline(projectdir("logs","progress",simulation)),"\r")[end][1:end-3]))
        end
        model.tableData[] = DataTable(df)
        sleep(1)
    end

    onany(model.simulate_button) do (_...)
        #Launch simulation
        println("Launching simulation")
        worker_path = srcdir("Dashboard","simulation_worker.jl")
        threads = TumorSim.Config["Worker"]["threads"]

        worker = `julia -p $threads --check-bounds=yes $worker_path 
                        $(model.death_rate[].range.start/1000) 0.05 $(model.death_rate[].range.stop/1000) 
                        $(model.mutation_rate[].range.start/1000) 0.005 $(model.mutation_rate[].range.stop/1000) 
                        10000 
                        $(join(model.scenario[],"")) 
                        10 
                        $(model.detecting_size[].range.start) 500 $(model.detecting_size[].range.stop)
                        $(model.starting_size[].range.start/1000) 0.05 $(model.starting_size[].range.stop/1000)
                        $(model.pausing_size[].range.start/1000) 0.05 $(model.pausing_size[].range.stop/1000)
                        $(model.kill_rate[].range.start/1000) 0.05 $(model.kill_rate[].range.stop/1000)
                        $(model.cost_of_resistance[].range.start/1000) 0.05 $(model.cost_of_resistance[].range.stop/1000)
                        $(model.repetitions[])
                        $(model.name[])`
        println(worker)
        #Aquí deberiamos de crear el archivo donde se va a guardar el log y poner algo como "Inicializando..." para que aparezca instantaneamente.
        #run(worker, wait=false)
        @async run(worker)
    end

    onany(model.cancel_button) do (_...)
        #Cancel simulation
        println("Canceling simulation")
    end

    onany(model.death_rate, model.mutation_rate, model.cost_of_resistance, model.detecting_size, model.starting_size, model.pausing_size, model.kill_rate, model.scenario, model.repetitions) do (_...)
        model.predicted_simulations[] = 2 * length(model.death_rate[].range.start:50:model.death_rate[].range.stop) * length(model.mutation_rate[].range.start:5:model.mutation_rate[].range.stop) *  length(model.cost_of_resistance[].range.start:50:model.cost_of_resistance[].range.stop) *  length(model.detecting_size[].range.start:500:model.detecting_size[].range.stop) *  length(model.starting_size[].range.start:50:model.starting_size[].range.stop) *  length(model.pausing_size[].range.start:50:model.pausing_size[].range.stop) *  length(model.kill_rate[].range.start:50:model.kill_rate[].range.stop) * length(model.scenario[]) * model.repetitions[]
    end

    page(model, class="container", title="TumorSim Simulation Dashboard",
        head_content=Genie.Assets.favicon_support(),
        prepend=style(
            """
            tr:nth-child(even) {
              background: #F8F8F8 !important;
            }
            .modebar {
              display: none!important;
            }
            .st-module {
              marign: 20px;
              background-color: #FFF;
              border-radius: 5px;
              box-shadow: 0px 4px 10px rgba(0, 0, 0, 0.04);
            }
            .stipple-core .st-module > h5,
            .stipple-core .st-module > h6 {
              border-bottom: 0px !important;
            }
            """
        ),
        [
            heading("TumorSim Simulation Dashboard")
            row([
                cell(
                class="st-module",
                [
                    btn("View results", color="primary", textcolor="black", type="a", href="view")
                    btn("Launch simulations", color="primary", textcolor="black", type="a", href="launch")
                ])
            ])
            
            row([
                cell(
                    class="st-module",
                    [
                    h5("Currently running simulations: ")
                    table(:tableData; pagination=:credit_data_pagination)
                ])
            ])
            row([
                cell(
                class="st-module",
                [
                    row([
                    btn("Launch a new simulation", color="primary", textcolor="black", @click("simulate_button += 1"))
                    h5(["  <---   Click to launch ",
                    GenieFramework.span(model.predicted_simulations, @text(:predicted_simulations)),
                    " simulations with the parameters below."])        
                    ])
                    br()
                    h5("Simulation name:  ")
                    input("", @bind(:name))
                    
                ]
                )
            ])

            row([
                cell(
                    class="st-module",
                    [
                        h6("Death rate (x1000)(Turnover)")
                        Stipple.range(0:50:1000,
                            @data(:death_rate);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Mutation rate (x1000)")
                        Stipple.range(0:5:200,
                            @data(:mutation_rate);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Cost of resistance (x1000)")
                        Stipple.range(0:50:1000,
                            @data(:cost_of_resistance);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Detecting size")
                        Stipple.range(0:500:5000,
                            @data(:detecting_size);
                            label=true)
                    ]
                )])
            row([
                cell(
                    class="st-module",
                    [
                        h6("Starting size % (x1000)")
                        Stipple.range(0:50:1000,
                            @data(:starting_size);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Pausing size % (x1000)")
                        Stipple.range(0:50:1000,
                            @data(:pausing_size);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Kill rate % (x1000)")
                        Stipple.range(0:50:1000,
                            @data(:kill_rate);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Scenario Dimensionality")
                        Stipple.select(:scenario; options=:scenario_list, multiple=true) 
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h5("Number of repetitions")
                        slider(1:1:100,
                            @data(:repetitions);
                            label=true)
                    ]
                    )
                ])
            

        ]
    )
end

launch_simulations_html = LaunchSimulationsPage |> init |> launch_simulations_ui |> html

route("/launch") do
    launch_simulations_html
end

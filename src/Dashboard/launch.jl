
using Stipple, StipplePlotly, StippleUI, Genie, GenieFramework
using DataFrames
using Statistics
using BSON


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

    pr::R{RangeData{Int64}} = RangeData(27:27)
    dr::R{RangeData{Int64}} = RangeData(500:500)
    mr::R{RangeData{Int64}} = RangeData(10:10)

    cr::R{RangeData{Int64}} = RangeData(200:200)

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

        worker = `julia -p 8 --check-bounds=yes $worker_path 
                        $(model.pr[].range.start/1000) 0.005 $(model.pr[].range.stop/1000) 
                        $(model.dr[].range.start/1000) 0.05 $(model.dr[].range.stop/1000) 
                        $(model.mr[].range.start/1000) 0.005 $(model.mr[].range.stop/1000) 
                        1000000 
                        $(join(model.scenario[],"")) 
                        10 
                        $(model.detecting_size[].range.start) 100 $(model.detecting_size[].range.stop)
                        $(model.starting_size[].range.start/1000) 0.05 $(model.starting_size[].range.stop/1000)
                        $(model.pausing_size[].range.start/1000) 0.05 $(model.pausing_size[].range.stop/1000)
                        $(model.kill_rate[].range.start/1000) 0.05 $(model.kill_rate[].range.stop/1000)
                        $(model.cr[].range.start/1000) 0.05 $(model.cr[].range.stop/1000)
                        $(model.repetitions[])
                        $(model.name[])`
        println(worker)
        @async run(worker)
        
    end

    onany(model.cancel_button) do (_...)
        #Cancel simulation
        println("Canceling simulation")
    end

    onany(model.dr, model.pr, model.mr, model.cr, model.detecting_size, model.starting_size, model.pausing_size, model.kill_rate, model.scenario, model.repetitions) do (_...)
        model.predicted_simulations[] = 2 * length(model.dr[].range.start:50:model.dr[].range.stop) * length(model.pr[].range.start:5:model.pr[].range.stop) * length(model.mr[].range.start:5:model.mr[].range.stop) *  length(model.cr[].range.start:50:model.cr[].range.stop) *  length(model.detecting_size[].range.start:100:model.detecting_size[].range.stop) *  length(model.starting_size[].range.start:50:model.starting_size[].range.stop) *  length(model.pausing_size[].range.start:50:model.pausing_size[].range.stop) *  length(model.kill_rate[].range.start:50:model.kill_rate[].range.stop) * length(model.scenario[]) * model.repetitions[]
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
                        h6("Proliferation rate (x1000)")
                        Stipple.range(0:10:200,
                            @data(:pr);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Death rate (x1000)(Turnover)")
                        Stipple.range(0:50:1000,
                            @data(:dr);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Mutation rate (x1000)")
                        Stipple.range(0:5:200,
                            @data(:mr);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Cost of resistance (x1000)")
                        Stipple.range(0:50:1000,
                            @data(:cr);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Detecting size")
                        Stipple.range(0:100:5000,
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

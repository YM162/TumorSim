
using Stipple, StipplePlotly, StippleUI, Genie, GenieFramework
using DataFrames
using Statistics
using BSON

@reactive mutable struct ViewResultsPage <: ReactiveModel

    button::R{Int} = 0
    clear_button::R{Int} = 0

    loading::R{Bool} = false
    disabled::R{Bool} = false
    
    visible_simulations::R{Int} = 0

    simulations::R{DataFrame} = DataFrame()

    scenario_list::R{Vector{Int64}} = [0, 1, 2, 3]
    scenario::R{Vector{Int64}} = [3]

    x_variable_list::R{Vector{Symbol}} = [:pr, :dr, :mr, :cr, :t_kill_rate, :t_detecting_size, :t_starting_size, :t_pausing_size, :s_dim]
    x_variable::R{Symbol} = :dr

    simulation_results_list::R{Vector{String}} = readdir(datadir("simulations"))
    simulation_results::R{String} = "Select a file to load"

    y_variable_list::R{Vector{Symbol}} = [:TTP, :recovery_rate]
    y_variable::R{Symbol} = :TTP
                                          
    detecting_size::R{RangeData{Int64}} = RangeData(0:5000)
    starting_size::R{RangeData{Int64}} = RangeData(0:1000)
    pausing_size::R{RangeData{Int64}} = RangeData(0:1000)

    pr::R{RangeData{Int64}} = RangeData(0:200)
    dr::R{RangeData{Int64}} = RangeData(0:1000)
    mr::R{RangeData{Int64}} = RangeData(0:200)

    cr::R{RangeData{Int64}} = RangeData(0:1000)

    kill_rate::R{RangeData{Int64}} = RangeData(0:1000)

    repetitions::R{Int64} = 10

    plot_data::R{Vector{PlotData}} = []
    layout::R{PlotLayout} = PlotLayout(
                                    plot_bgcolor = "#FFFFFF",
                                    title = PlotLayoutTitle(text="Simulation results", font=Font(24)),
                                    xaxis = [PlotLayoutAxis(xy="x",title_text="dr")],
                                    yaxis = [PlotLayoutAxis(xy="y",title_text="TTP")],
                                    )

    config::R{PlotConfig} = PlotConfig()

end


function tumorPlot(model::ViewResultsPage,name) 
    if name == "Continuous therapy"
        data = filter((row) -> row["t_pausing_size"] == 0.0,model.simulations[])

    elseif name == "Adaptive therapy"
        data = filter((row) -> row["t_pausing_size"] != 0,model.simulations[]) #Hay que añadir la t_pausing_size a las variables para cambiarlo dependiendo del detecting size
    end

    #We apply all Filters

    filter!((row) -> (row["pr"] >= (model.pr[].range.start/1000)) & (row["pr"] <= (model.pr[].range.stop/1000)), data)
    filter!((row) -> (row["dr"] >= (model.dr[].range.start/1000)) & (row["dr"] <= (model.dr[].range.stop/1000)), data)
    filter!((row) -> (row["mr"] >= (model.mr[].range.start/1000)) & (row["mr"] <= (model.mr[].range.stop/1000)), data)
    filter!((row) -> (row["cr"] >= (model.cr[].range.start/1000)) & (row["cr"] <= (model.cr[].range.stop/1000)), data)
    filter!((row) -> (row["t_kill_rate"] >= (model.kill_rate[].range.start/1000)) & (row["t_kill_rate"] <= (model.kill_rate[].range.stop/1000)), data)
    filter!((row) -> (row["t_starting_size"] >= (model.starting_size[].range.start/1000)) & (row["t_starting_size"] <= (model.starting_size[].range.stop/1000)), data)
    filter!((row) -> (row["t_detecting_size"] >= (model.detecting_size[].range.start)) & (row["t_detecting_size"] <= (model.detecting_size[].range.stop)), data)
    filter!((row) -> (row["t_pausing_size"] >= (model.pausing_size[].range.start/1000)) & (row["t_pausing_size"] <= (model.pausing_size[].range.stop/1000)), data)
    filter!((row) -> (row["s_dim"] in model.scenario[]), data)

    model.visible_simulations[] = length(data[!,"TTP"]) * 2

    #We continue with the plot
    x = sort(unique(data[!, string(model.x_variable[])]))
    y = []
    for x_ in x
        x_filter = filter((row) -> row[string(model.x_variable[])] == x_ ,data)
        if model.y_variable[] == :TTP
            y_ = filter((row) -> row[string(:TTP)] != -1 ,x_filter)
            if length(y_[!,"TTP"]) == 0
                append!(y,-1)
            else
                append!(y,mean(y_[!, string(:TTP)]))
            end
        elseif model.y_variable[] == :recovery_rate
            y_ = filter((row) -> row[string(:TTP)] == -1 ,x_filter)
            append!(y,length(y_[!,"TTP"])/length(x_filter[!,"TTP"]))
        end
    end

    model.layout[] = PlotLayout(
                                plot_bgcolor = "#FFFFFF",
                                title = PlotLayoutTitle(text="Simulation results", font=Font(24)),
                                xaxis = [PlotLayoutAxis(xy="x",title_text=string(model.x_variable[]))],
                                yaxis = [PlotLayoutAxis(xy="y",title_text=string(model.y_variable[]))],
                                )

    PlotData(
        x=x,
        y = y,
        plot=StipplePlotly.Charts.PLOT_TYPE_SCATTER,
        name=name,
    )
end




function ui(model::ViewResultsPage)

    onany(model.button) do (_...)
        #Load data
        if model.simulation_results[] == "Select a file to load"
            return
        end

        model.loading[] = true
        filepath = datadir("simulations", model.simulation_results[])
 
        append!(model.simulations[],BSON.load(filepath)["df"])

        unique!(model.simulations[])
        
        model.plot_data[] = [tumorPlot(model,"Continuous therapy"),tumorPlot(model,"Adaptive therapy")]
        model.loading[] = false
    end

    onany(model.clear_button) do (_...)
        model.simulations[] = DataFrame()
        model.visible_simulations[] = 0
        model.plot_data[] = []
    end


    onany(model.x_variable, model.y_variable, model.pr, model.dr, model.mr, model.cr, model.detecting_size, model.starting_size, model.pausing_size, model.kill_rate,model.scenario) do (_...)
        model.plot_data[] = [tumorPlot(model,"Continuous therapy"),tumorPlot(model,"Adaptive therapy")]
    end

    onany(model.simulations) do (_...)
        unique!(model.simulations[])
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
                btn("Load data", loading=@data(:loading),disabled=@data(:disabled),color="primary", textcolor="black", @click("button += 1"), [
                    tooltip(contentclass="bg-indigo", contentstyle="font-size: 16px",
                        style="offset: 10px 10px", "Load all simulations")])
                cell(
                    class="st-module",
                    [
                        Stipple.select(:simulation_results; options=:simulation_results_list)
                    ])
                btn("Clear ALL data",color="primary", textcolor="black", @click("clear_button += 1"), [
                    tooltip(contentclass="bg-indigo", contentstyle="font-size: 16px",
                        style="offset: 10px 10px", "Load all simulations")])
            ])
            
            row([
                cell(
                    size=6,
                    class="st-module",
                    [   
                        row([
                            cell(
                                class="st-module",
                                [
                                    Stipple.select(:y_variable; options=:y_variable_list)
                                ])
                            cell(
                                class="st-module",
                                [
                                    h5("by")
                                ])
                            cell(
                                class="st-module",
                                [
                                    Stipple.select(:x_variable; options=:x_variable_list)
                                ])
                        ])
                        
                        plot(:plot_data, layout=:layout, config=:config)
                        
                        ]
                )
                cell(
                    class="st-module",
                    [
                        cell(
                            class="st-module",
                            [
                                h5(["Filters: (Currently showing ",
                                GenieFramework.span(model.visible_simulations, @text(:visible_simulations)),
                                " unique simulations.)"])
                                br()
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
                        )
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
                    ]
                    
                )
            ])
        ]
    )
end

htmlfile = ViewResultsPage |> init |> ui |> html

route("/view") do
    htmlfile
end

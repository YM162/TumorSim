#Build with a lot of references from https://github.com/GenieFramework/StippleDemos

using Distributed

@everywhere using Stipple, StipplePlotly, StippleUI
@everywhere using DataFrames
@everywhere using Statistics

@reactive mutable struct HeatPages <: ReactiveModel


    simulate_text::R{String} = "Simulate!"
    value::R{Int} = 0
    click::R{Int} = 0

    completed_simulations::R{Int} = 0
    total_simulations::R{Int} = 4320

    scenario_list::R{Vector{Symbol}} = [:_0D, :_1D, :_2D, :_3D]
    scenario::R{Symbol} = :_2D

    x_variable_list::R{Vector{Symbol}} = [:pr, :dr, :mr]
    x_variable::R{Symbol} = :dr

    y_variable_list::R{Vector{Symbol}} = [:TTP, :recovery_rate]
    y_variable::R{Symbol} = :TTP

    fitness::R{Dict{Vector{Int64},Real}} = Dict([0,0,0]=>1, 
                                                [1,0,0]=>1.3,
                                                [0,1,0]=>1.2,
                                                [1,1,0]=>1.5,
                                                [1,1,1]=>1.2)

    df::R{DataFrame} = DataFrame(mr=Float64[], pr=Float64[], dr=Float64[], TTP=Int[], t_pausing_size=Int[])
                                                
                                
    pr::R{RangeData{Int64}} = RangeData(20:35)
    dr::R{RangeData{Int64}} = RangeData(350:750)
    mr::R{RangeData{Int64}} = RangeData(0:25)

    repetitions::R{Int64} = 10

    plot_data::R{Vector{PlotData}} = []
    layout::R{PlotLayout} = PlotLayout(
                                    plot_bgcolor = "#FFFFFF",
                                    title = PlotLayoutTitle(text="TTP", font=Font(24))
                                    )

    config::R{PlotConfig} = PlotConfig()

end


function tumorPlot(model::HeatPages,name) 
    if name == "Continuous therapy"
        data = filter((row) -> row["t_pausing_size"] ==0,model.df[])
    elseif name == "Adaptive therapy"
        data = filter((row) -> row["t_pausing_size"] ==1000,model.df[]) #Hay que añadir la t_pausing_size a las variables para cambiarlo dependiendo del detecting size
    end

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

    PlotData(
        x=x,
        y = y,
        plot=StipplePlotly.Charts.PLOT_TYPE_SCATTER,
        name=name,
    )
end




function ui(model::HeatPages)

    onany(model.value) do (_...)
        model.simulate_text[] = "Simulating..."
        model.click[] += 1

        scenario = model.scenario[]
        if scenario == :_0D
            scenario = create_scenario((1000000),10)
        elseif scenario == :_1D
            scenario = create_scenario((1000000,),10,"center")
        elseif scenario == :_2D
            scenario = create_scenario((1000,1000),10,"center")
        elseif scenario == :_3D
            scenario = create_scenario((100,100,100),10,"center")
        end

        parameters = Dict(
            "pr" => collect(model.pr[].range.start:5:model.pr[].range.stop)/1000,
            "dr" => collect(model.dr[].range.start:50:model.dr[].range.stop)/1000,
            "mr" => collect(model.mr[].range.start:5:model.mr[].range.stop)/1000,   
            "scenario" => scenario, 
            "fitness" => model.fitness[],
            "treatment" => [create_treatment(3000, 2000, 1000, 3, 0.75) ,create_treatment(3000, 2000, 0, 3, 0.75) ],
            "seed" => map(abs,rand(Int64,model.repetitions[]))
        )

        parameter_combinations = dict_list(parameters)

        
        
        println("Starting simulation with $(length(parameter_combinations)) parameter combinations")
        model.completed_simulations[] = 0

        for chunk in Iterators.partition(parameter_combinations, 10)
            result_list = pmap(simulate,chunk,fill(3000,length(chunk)))
            for results in result_list
                push!(model.df[], [results["mr"], results["pr"], results["dr"], results["TTP"], results["t_pausing_size"]])
            end
            model.completed_simulations[] += length(chunk)
            model.plot_data[] = [tumorPlot(model,"Continuous therapy"),tumorPlot(model,"Adaptive therapy")]
        end

        model.plot_data[] = [tumorPlot(model,"Continuous therapy"),tumorPlot(model,"Adaptive therapy")]
        println(model.df[])
        model.simulate_text[] = "Simulate!"
    end

    onany(model.dr, model.pr, model.mr, model.scenario, model.repetitions) do (_...)
        model.total_simulations[] = 2 * length(model.dr[].range.start:50:model.dr[].range.stop) * length(model.pr[].range.start:5:model.pr[].range.stop) * length(model.mr[].range.start:5:model.mr[].range.stop) * model.repetitions[]
    end

    onany(model.x_variable, model.y_variable) do (_...)
        model.plot_data[] = [tumorPlot(model,"Continuous therapy"),tumorPlot(model,"Adaptive therapy")]
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
                        h6("Proliferation rate (x1000)")
                        Stipple.range(10:5:100,
                            @data(:pr);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Death rate (x1000)")
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
                        h6("Repetitions")
                        slider(1:1:100,
                            @data(:repetitions);
                            label=true)
                    ]
                )
                cell(
                    class="st-module",
                    [
                        h6("Scenario Dimensionality")
                        Stipple.select(:scenario; options=:scenario_list)
                    ]
                )])
            row([
                btn(model.simulate_text[], color="primary", textcolor="black", @click("value += 1"), [
                    tooltip(contentclass="bg-indigo", contentstyle="font-size: 16px",
                        style="offset: 10px 10px", "Click the button to start simulation")])
                cell(
                    class="st-module",
                    [
                        h6(["Simulation progress: ",
                            GenieFramework.span(model.completed_simulations, @text(:completed_simulations)),
                            "/",
                            GenieFramework.span(model.total_simulations, @text(:total_simulations))
                            ])
                    ])
            ])
            row([
                cell(
                    size=6,
                    class="st-module",
                    [
                        h5("Result Plot")
                        plot(:plot_data, layout=:layout, config=:config)
                        
                        ]
                )
                cell(
                    class="st-module",
                    [
                        h5("X Variable")
                        Stipple.select(:x_variable; options=:x_variable_list)
                        h5("Y Variable")
                        Stipple.select(:y_variable; options=:y_variable_list)
                    ]
                    
                )
            ])
        ]
    )
end

htmlfile = HeatPages |> init |> ui |> html

route("/") do
    htmlfile
end

up()
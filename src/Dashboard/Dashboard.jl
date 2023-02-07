#Todo esto hay que meterlo en un submodulo como en todos los demás.
module Dashboard
    export launch_dashboard, kill_dashboard
    using DrWatson

    using FilePathsBase
    using Stipple, StipplePlotly, StippleUI, Genie, GenieFramework
    using DataFrames
    using Statistics
    using BSON

    using TumorSim.Fitness
    using TumorSim.Analysis
    using TumorSim.Scenario
    using TumorSim.Simulate
    using TumorSim.Treatment
    using TumorSim.TumorModel
    using TumorSim
    function launch_dashboard()

        include(srcdir("Dashboard", "view.jl"))
        include(srcdir("Dashboard", "launch.jl"))

        route("/") do
            Genie.Renderer.redirect("/view")
        end

        up(TumorSim.Config["Dashboard"]["port"],TumorSim.Config["Dashboard"]["ip"])
    end

    function kill_dashboard()
        Genie.down!()
    end
end
using .Dashboard
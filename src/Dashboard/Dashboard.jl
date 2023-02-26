#Todo esto hay que meterlo en un submodulo como en todos los demás.
module Dashboard
    export launch_dashboard, kill_dashboard
    using DrWatson

    using FilePathsBase
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

    end

end
using .Dashboard
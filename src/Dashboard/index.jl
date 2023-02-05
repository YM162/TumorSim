#Built with a lot of references from https://github.com/GenieFramework/StippleDemos

using Stipple, StipplePlotly, StippleUI, Genie, GenieFramework

include(srcdir("Dashboard", "view.jl"))
#include(srcdir("Dashboard", "launch.jl"))

route("/") do
    Genie.Renderer.redirect("/view")
end

up()
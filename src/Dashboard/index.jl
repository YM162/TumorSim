using GenieFramework
@genietools

@handlers begin
    @out message = "Number of characters:"
    @in text = "Write here"

    @onchange isready begin
        @show "App is loaded"
    end
    
    @onchange text begin
        message = "Number of characters: $(length(text))"
    end
end

@page("/", srcdir("Dashboard/index.jl.html"))

up()
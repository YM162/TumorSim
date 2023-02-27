
#We define the scenario and the functions to create it using multiple dispatch.
module Scenario
    export ScenarioObject, create_scenario
    struct ScenarioObject
        x::Int64
        y::Int64
        z::Int64
        cell_pos::Vector{Tuple}
        mix::Bool
    end
    
    #We return the same scenario but we can only have one cell in each position. We move each of the cells to the nearest empty space without going out of the space limits.
    function flatten_scenario(scenario::ScenarioObject)
        restricted_cells = []
        dimensions = [scenario.x,scenario.y,scenario.z]
        for cell in scenario.cell_pos
            newpos = collect(cell)
            timeout = 0
            while Tuple(newpos) in restricted_cells
                for dim in eachindex(newpos)
                    newpos[dim] = newpos[dim]+rand([-1,0,1])
                    if newpos[dim] > dimensions[dim] || newpos[dim] < 1
                        newpos[dim] = dimensions[dim]
                    end
                end
                if timeout > 1000
                    print("Error, there is not enough space for every cell.")
                    return
                end
                timeout = timeout + 1
            end
            newpos = Tuple(newpos)
            push!(restricted_cells,newpos)
        end
        return ScenarioObject(scenario.x,scenario.y,scenario.z,restricted_cells,scenario.mix)
    end

    #1D
    function create_scenario(size::Tuple{Int64},ncells::Int,cell_pos::String="center",mix::Bool=false)
        if cell_pos=="random"
            cell_pos = [((rand(1:size[1])),1,1) for _ in 1:ncells]
        elseif cell_pos=="center"
            cell_pos = [(Int(floor(size[1]/2))+1,1,1) for _ in 1:ncells]
        else
            print("Error, cell_pos option not found.")
            return
        end
            return flatten_scenario(ScenarioObject(size[1],1,1,cell_pos,mix))
    end

    function create_scenario(size::Tuple{Int64},cell_pos::Vector{Tuple{Int64}},mix::Bool=false)
            return flatten_scenario(ScenarioObject(size[1],1,1,[(x[1],1,1) for x in cell_pos],mix))
    end

    #2D
    function create_scenario(size::Tuple{Int64, Int64},ncells::Int,cell_pos::String="center",mix::Bool=false)
        if cell_pos=="random"
            cell_pos = [((rand(1:size[1])),(rand(1:size[2])),1) for _ in 1:ncells]
        elseif cell_pos=="center"
            cell_pos = [(Int(floor(size[1]/2))+1,Int(floor(size[2]/2))+1,1) for _ in 1:ncells]
        else
            print("Error, cell_pos option not found.")
            return
        end
            return flatten_scenario(ScenarioObject(size[1],size[2],1,cell_pos,mix))
    end

    function create_scenario(size::Tuple{Int64,Int64},cell_pos::Vector{Tuple{Int64, Int64}},mix::Bool=false)
            return flatten_scenario(ScenarioObject(size[1],size[2],1,[(x[1],x[2],1) for x in cell_pos],mix))
    end

    #3D
    function create_scenario(size::Tuple{Int64, Int64, Int64},ncells::Int,cell_pos::String="center",mix::Bool=false)
        if cell_pos=="random"
            cell_pos = [((rand(1:size[1])),(rand(1:size[2])),(rand(1:size[3]))) for _ in 1:ncells]
        elseif cell_pos=="center"
            cell_pos = [(Int(floor(size[1]/2))+1,Int(floor(size[2]/2))+1,Int(floor(size[3]/2))+1) for _ in 1:ncells]
        else
            print("Error, cell_pos option not found.")
            return
        end
            return flatten_scenario(ScenarioObject(size[1],size[2],size[3],cell_pos,mix))
    end
    π
    function create_scenario(size::Tuple{Int64,Int64,Int64},cell_pos::Vector{Tuple{Int64, Int64,Int64}},mix::Bool=false)
            return flatten_scenario(ScenarioObject(size[1],size[2],size[3],[(x[1],x[2],x[3]) for x in cell_pos],mix))
    end
end
using .Scenario
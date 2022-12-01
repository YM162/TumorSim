
#We define the scenario and the functions to create it using multiple dispatch.
struct Scenario
    x::Int
    y::Int
    z::Int
    cell_pos::Vector
    wall_pos::Vector
end

#0D
function create_scenario(size::Int64,ncells::Int)   
    return Scenario(size,0,0,[(1,1,1) for i in 1:ncells],[])
end

#1D
function create_scenario(size::Tuple{Int64},ncells::Int,cell_pos::String="center",wall_pos::Vector=[])
    if cell_pos=="random"
        cell_pos = [((rand(1:size[1])),1,1) for i in 1:ncells]
    elseif cell_pos=="center"
        cell_pos = [(Int(floor(size[1]/2))+1,1,1) for i in 1:ncells]
    else
        print("Error, cell_pos option not found.")
        return
    end
        return Scenario(size[1],1,1,cell_pos,wall_pos)
end

function create_scenario(size::Tuple{Int64},cell_pos::Vector{Tuple{Int64}},wall_pos=[])
        return Scenario(size[1],1,1,[(x[1],1,1) for x in cell_pos],wall_pos)
end

#2D
function create_scenario(size::Tuple{Int64, Int64},ncells::Int,cell_pos::String="center",wall_pos::Vector=[])
    if cell_pos=="random"
        cell_pos = [((rand(1:size[1])),(rand(1:size[2])),1) for i in 1:ncells]
    elseif cell_pos=="center"
        cell_pos = [(Int(floor(size[1]/2))+1,Int(floor(size[2]/2))+1,1) for i in 1:ncells]
    else
        print("Error, cell_pos option not found.")
        return
    end
        return Scenario(size[1],size[2],1,cell_pos,wall_pos)
end

function create_scenario(size::Tuple{Int64,Int64},cell_pos::Vector{Tuple{Int64, Int64}},wall_pos::Vector=[])
        return Scenario(size[1],size[2],1,[(x[1],x[2],1) for x in cell_pos],wall_pos)
end

#3D
function create_scenario(size::Tuple{Int64, Int64, Int64},ncells::Int,cell_pos::String="center",wall_pos::Vector=[])
    if cell_pos=="random"
        cell_pos = [((rand(1:size[1])),(rand(1:size[2])),(rand(1:size[3]))) for i in 1:ncells]
    elseif cell_pos=="center"
        cell_pos = [(Int(floor(size[1]/2))+1,Int(floor(size[2]/2))+1,Int(floor(size[3]/2))+1) for i in 1:ncells]
    else
        print("Error, cell_pos option not found.")
        return
    end
        return Scenario(size[1],size[2],size[3],cell_pos,wall_pos)
end

function create_scenario(size::Tuple{Int64,Int64,Int64},cell_pos::Vector{Tuple{Int64, Int64,Int64}},wall_pos::Vector=[])
        return Scenario(size[1],size[2],size[3],[(x[1],x[2],x[3]) for x in cell_pos],wall_pos)
end
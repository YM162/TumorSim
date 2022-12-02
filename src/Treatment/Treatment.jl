#We define what a treatment is
mutable struct Treatment
    detecting_size::Int
    starting_size::Int
    pausing_size::Int 
    resistance_gene::Int
    kill_rate::Float64
    active::Bool
    detected::Bool
end

function create_treatment(detecting_size,starting_size,pausing_size,resistance_gene,kill_rate)
    return Treatment(detecting_size,starting_size,pausing_size,resistance_gene,kill_rate,false,false)
end
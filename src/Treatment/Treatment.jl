#We define what a treatment is
mutable struct Treatment
    detecting_size::Int
    starting_size::Int
    pausing_size::Int 
    resistance_gene::Int
    kill_rate::Float16
    active::Bool
    detected::Bool
end


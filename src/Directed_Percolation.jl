#=
Directed_Percolation.jl

This module provides the core functions for simulating a (1+1)-dimensional 
directed percolation process. It is designed to be imported into other scripts
or notebooks for analysis and visualization.

Core Logic:
- `evolve`: The main function to run the time evolution of the system.
- `generate_initial_state`: A helper to create random initial conditions.
- `calculate_average_density`: A helper to compute the order parameter.
=#

module DirectedPercolation

using Random
using Statistics

# Export the functions that will be used by other programs
export evolve, generate_initial_state, calculate_average_density

"""
    evolve(N::Int, p::Float64, q::Float64, t_max::Int, initial_state::Vector{Int})

Evolves the directed percolation system for `t_max` steps.

# Arguments
- `N`: The size of the lattice (number of sites).
- `p`: The base probability parameter.
- `q`: The survival rate parameter.
- `t_max`: The total number of time steps to simulate.
- `initial_state`: A vector of 0s and 1s representing the state at t=0.

# Returns
- A `t_max+1` x `N` matrix where each row is the state of the system at a given time step.
"""
function evolve(N::Int, p::Float64, q::Float64, t_max::Int, initial_state::Vector{Int})
    # Check for valid inputs
    if !(0.0 <= p <= 1.0) || !(0.0 <= q <= 1.0)
        error("Probabilities p and q must be between 0 and 1.")
    end
    if length(initial_state) != N
        error("Length of initial_state must be equal to N.")
    end

    # The effective probability after considering survival rate
    pq = p * q

    # Probabilities for different neighbor configurations
    prob_one_neighbor = pq
    prob_two_neighbors = 2 * pq - pq^2

    # History matrix to store the state at each time step
    history = zeros(Int, t_max + 1, N)
    history[1, :] = initial_state

    # Time evolution loop
    for t in 1:t_max
        current_state = @view history[t, :]
        next_state = @view history[t+1, :]
        for i in 1:N
            # Periodic boundary conditions: left neighbor is at i-1, or N if i is 1.
            left_neighbor_idx = (i == 1) ? N : i - 1
            
            n_left = current_state[left_neighbor_idx]
            n_current = current_state[i]
            
            num_active_neighbors = n_left + n_current

            prob_active = 0.0
            if num_active_neighbors == 1
                prob_active = prob_one_neighbor
            elseif num_active_neighbors == 2
                prob_active = prob_two_neighbors
            end

            # Update the site based on the calculated probability
            if rand() < prob_active
                next_state[i] = 1
            else
                next_state[i] = 0
            end
        end
    end
    return history
end

"""
    generate_initial_state(N::Int; density::Float64=0.5)

Generates a random initial state with a given density of active sites.
"""
function generate_initial_state(N::Int; density::Float64=0.5)
    return rand(N) .< density
end

"""
    calculate_average_density(history::Matrix{Int})

Calculates the average density of active sites for each time step.
"""
function calculate_average_density(history::Matrix{Int})
    N = size(history, 2)
    return vec(mean(history, dims=2))
end

end # end of module DirectedPercolation

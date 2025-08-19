#=
State_Evol.jl

This file contains the core functions for evolving the state of the 
directed percolation system, implementing a staggered time evolution rule.
=#

"""
    evolve(N::Int, p::Float64, q::Float64, t_max::Int, initial_state::Vector{Int})

Evolves the directed percolation system for `t_max` steps using a staggered
update rule that depends on whether the time step is even or odd.

- On even time steps `t`, site `i` at `t+1` is influenced by sites `i-1` and `i` at `t`.
- On odd time steps `t`, site `i` at `t+1` is influenced by sites `i` and `i+1` at `t`.
"""
function evolve(N::Int, p::Float64, q::Float64, t_max::Int, initial_state::Vector{Int})
    # Input validation
    if !(0.0 <= p <= 1.0) || !(0.0 <= q <= 1.0)
        error("Probabilities p and q must be between 0 and 1.")
    end
    if length(initial_state) != N
        error("Length of initial_state must be equal to N.")
    end

    pq = p * q
    prob_one_neighbor = pq
    prob_two_neighbors = 2 * pq - pq^2

    history = zeros(Int, t_max + 1, N)
    history[1, :] = initial_state

    # --- Main Time Evolution Loop ---
    for t in 1:t_max
        current_state = @view history[t, :]
        next_state = @view history[t+1, :]

        # Loop over each site in the lattice
        for i in 1:N
            
            # --- Determine neighbors based on even/odd time step ---
            local n1, n2 # Neighbor states
            if t % 2 == 0 # Current time step t is EVEN
                # Next state at site `i` depends on current state at `i-1` and `i`
                left_neighbor_idx = (i == 1) ? N : i - 1
                n1 = current_state[left_neighbor_idx]
                n2 = current_state[i]
            else # Current time step t is ODD
                # Next state at site `i` depends on current state at `i` and `i+1`
                right_neighbor_idx = (i == N) ? 1 : i + 1
                n1 = current_state[i]
                n2 = current_state[right_neighbor_idx]
            end

            num_active_neighbors = n1 + n2

            prob_active = 0.0
            if num_active_neighbors == 1
                prob_active = prob_one_neighbor
            elseif num_active_neighbors == 2
                prob_active = prob_two_neighbors
            end

            # Update the site for the next time step
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
    return Vector{Int}(rand(N) .< density)
end

"""
    calculate_average_density(history::Matrix{Int})

Calculates the average density of active sites for each time step.
"""
function calculate_average_density(history::Matrix{Int})
    N = size(history, 2)
    return vec(mean(history, dims=2))
end

#=
Directed_Percolation.jl

This is the main module file for the DirectedPercolation package.
It loads dependencies, includes the core logic and plotting functions
from other files, and exports the public API.
=#

module DirectedPercolation

# --- Dependencies ---
using Random
using Statistics
using Plots

# --- Export Public API ---
# Core functions
export evolve, generate_initial_state, calculate_average_density
# Plotting functions
export plot_density_vs_time, plot_phase_diagram, plot_lattice_evolution

# --- Include source files ---
include("State_Evol.jl")
include("Plots_Helpers.jl")

end # end of module DirectedPercolation

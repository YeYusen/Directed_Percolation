#=
Plots_Helpers.jl

This file contains all helper functions for plotting and visualizing
the results of the directed percolation simulation.
=#

"""
    plot_density_vs_time(N::Int, p::Float64, q::Float64, t_max::Int, num_trials::Int)

Runs multiple simulations and plots the average active site density vs. time.
"""
function plot_density_vs_time(N::Int, p::Float64, q::Float64, t_max::Int, num_trials::Int)
    println("Generating plot for p=$p, q=$q...")
    total_density_over_time = zeros(Float64, t_max + 1)

    for _ in 1:num_trials
        initial_state = generate_initial_state(N, density=0.5)
        history = evolve(N, p, q, t_max, initial_state)
        density_over_time = calculate_average_density(history)
        total_density_over_time .+= density_over_time
    end

    avg_density = total_density_over_time / num_trials

    plot(0:t_max, avg_density,
        title="Average Active Site Density vs. Time\n(p=$p, q=$q, N=$N, trials=$num_trials)",
        xlabel="Time (t)",
        ylabel="Average Density",
        label="Density",
        legend=:topright,
        linewidth=2,
        dpi=300
    )
end

"""
    plot_phase_diagram(N::Int, t_final::Int, p_steps::Int, q_steps::Int, num_trials::Int)

Generates a 2D color-coded phase diagram of final active site density.
"""
function plot_phase_diagram(N::Int, t_final::Int, p_steps::Int, q_steps::Int, num_trials::Int)
    println("Generating phase diagram...")
    p_range = range(0, 1, length=p_steps)
    q_range = range(0, 1, length=q_steps)
    
    final_densities = zeros(Float64, length(q_range), length(p_range))

    progress_step = ceil(Int, p_steps / 10)

    for (j, p) in enumerate(p_range)
        for (i, q) in enumerate(q_range)
            avg_final_density = 0.0
            for _ in 1:num_trials
                initial_state = generate_initial_state(N, density=0.8)
                history = evolve(N, p, q, t_final, initial_state)
                density_over_time = calculate_average_density(history)
                start_index = floor(Int, 0.9 * t_final) + 1
                avg_final_density += mean(density_over_time[start_index:end])
            end
            final_densities[i, j] = avg_final_density / num_trials
        end
        
        if j % progress_step == 0 || j == p_steps
            println("Progress: $(round(j/p_steps*100, digits=1))%")
        end
    end

    heatmap(p_range, q_range, final_densities,
        title="Directed Percolation Phase Diagram\n(N=$N, t_final=$t_final, trials=$num_trials)",
        xlabel="p (Probability)",
        ylabel="q (Survival Rate)",
        c=:viridis,
        colorbar_title="Final Active Site Density",
        dpi=300,
        top_margin=5Plots.mm
    )
end

"""
    plot_staggered_lattice_grid(N::Int, t_max::Int; dpi=300, yflip_ax::Bool=true)

Draws the underlying staggered bond lattice efficiently by batching draw calls.
"""
function plot_staggered_lattice_grid(N::Int, t_max::Int; dpi=300, yflip_ax::Bool=true)
    
    # --- OPTIMIZATION: Collect all line coordinates before plotting ---
    # We use NaN to separate line segments in a single plot call.
    lines_x = Float64[]
    lines_y = Float64[]
    
    # Pre-allocate memory for slight performance gain
    sizehint!(lines_x, N * t_max * 2 * 3)
    sizehint!(lines_y, N * t_max * 2 * 3)

    for t in 0:t_max-1
        for i in 1:N
            x_start = (t % 2 == 0) ? i : i + 0.5
            y_start = Float64(t)
            y_end = Float64(t + 1)
            
            # Bond 1 (to the left)
            x_end1 = x_start - 0.5
            push!(lines_x, x_start, x_end1, NaN)
            push!(lines_y, y_start, y_end, NaN)

            # Bond 2 (to the right)
            x_end2 = x_start + 0.5
            push!(lines_x, x_start, x_end2, NaN)
            push!(lines_y, y_start, y_end, NaN)
        end
    end
    
    # --- Create the plot ---
    p = plot(
        title="Staggered Lattice Structure (N=$N, t_max=$t_max)",
        xlabel="Site",
        ylabel="Time",
        legend=false,
        grid=false,
        background_color=:black,
        dpi=dpi
    )

    # --- Single, efficient plot call for all bonds ---
    plot!(p, lines_x, lines_y, color=:gray40, lw=0.5)

    if yflip_ax
        yflip!(p)
    end
    
    return p
end


"""
    plot_lattice_evolution(history::Matrix{Int}; xlims=nothing, ylims=nothing, dpi=300, yflip_ax::Bool=true)

Generates a 2D scatter plot visualizing the state of the lattice over time.
- Staggers the lattice sites by 0.5 on odd time steps.
- `yflip_ax`: If true, inverts the y-axis so t=0 is at the top.
"""
function plot_lattice_evolution(history::Matrix{Int}; xlims=nothing, ylims=nothing, dpi=300, yflip_ax::Bool=true)
    t_max = size(history, 1) - 1
    N = size(history, 2)

    active_indices = findall(history .== 1)
    inactive_indices = findall(history .== 0)
    
    # Stagger the lattice: add 0.5 to x-coordinate for odd time steps
    stagger(site, time) = (time % 2 == 1) ? (site + 0.5, time) : (Float64(site), Float64(time))
    
    active_sites = [stagger(idx[2], idx[1] - 1) for idx in active_indices]
    inactive_sites = [stagger(idx[2], idx[1] - 1) for idx in inactive_indices]

    final_xlims = xlims === nothing ? (1, N) : xlims
    final_ylims = ylims === nothing ? (0, t_max) : ylims

    plot(
        title="Lattice Evolution (N=$N, t_max=$t_max)",
        xlabel="Site",
        ylabel="Time",
        legend=false,
        grid=false,
        background_color=:black,
        xlims=final_xlims,
        ylims=final_ylims,
        dpi=dpi
    )

    # Plot inactive sites with better visibility
    scatter!(first.(inactive_sites), last.(inactive_sites),
        marker=:circle,
        markercolor=:gray30, # Increased visibility
        markersize=2.0,     # Increased size
        markerstrokewidth=0
    )

    # Plot active sites
    scatter!(first.(active_sites), last.(active_sites),
        marker=:circle,
        markercolor=:cyan,
        markersize=2.0,     # Increased size
        markerstrokewidth=0
    )

    if yflip_ax
        yflip!(true)
    end
end

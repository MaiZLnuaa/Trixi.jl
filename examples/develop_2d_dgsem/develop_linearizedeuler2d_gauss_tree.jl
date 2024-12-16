#=
 *  @File          develop_linearizedeuler2d_gauss_tree.jl
 *
 *  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
 *  @Date          Mon Dec 16 2024 11:25:22 AM
 *
 *  @Description   1
 =#

using OrdinaryDiffEq
using Trixi

###############################################################################
# semidiscretization of the linearized Euler equations

equations = LinearizedEulerEquations2D(v_mean_global = (0.5, 0.0),c_mean_global = 1.0, rho_mean_global = 1.0)

# Create DG solver with polynomial degree = 5 and upwind flux as surface flux
solver = DGSEM(polydeg = 5, surface_flux = flux_godunov)

coordinates_min = (-100.0, -100.0) # minimum coordinates (min(x), min(y))
coordinates_max = (100.0, 100.0) # maximum coordinates (max(x), max(y))

# Create a uniformly refined mesh with periodic boundaries
mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level = 6, n_cells_max = 100_000, periodicity = false)

function initial_condition_gauss_wall(x, t, equations::LinearizedEulerEquations2D)
    v1_prime = 0.0
    v2_prime = 0.0
    rho_prime = p_prime = exp(-log(2)/log(ℯ)/9.0*(x[1]^2+x[2]^2))
    return SVector(rho_prime, v1_prime, v2_prime, p_prime)
end
initial_condition = initial_condition_gauss_wall

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions = boundary_condition_wall)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to 30.0
tspan = (0.0, 120.0)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval = 100)

# The SaveRestartCallback allows to save a file from which a Trixi.jl simulation can be restarted
# save_restart = SaveRestartCallback(interval = 100, save_final_restart = true)

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval = 100)
# save_solution = SaveSolutionCallback(dt = 10)

# The TimeSeriesCallback records the solution at the given points over time
time_series = TimeSeriesCallback(semi, [(0.0, 0.0), (-1.0, 0.5)])

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 1.0)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution, stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
sol = solve(ode, CarpenterKennedy2N54(thread = OrdinaryDiffEq.True(), williamson_condition = false),
dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks)

# Print the timer summary
summary_callback()
#=
 *  @File          develop_linearizedeuler_gauss_tree.jl
 *
 *  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
 *  @Date          Tue Dec 10 2024 3:16:28 PM
 *
 *  @Description   3d Gaussian pause
 =#

 using OrdinaryDiffEq
 using Trixi

###############################################################################
# semidiscretization of the linearized Euler equations

equations = LinearizedEulerEquations3D(v_mean_global = (0.0, 0.0, 0.0), c_mean_global = 1.0,
rho_mean_global = 1.0)

solver = DGSEM(polydeg = 5, surface_flux = flux_lax_friedrichs)

coordinates_min = (-5.0, -5.0, -5.0)
coordinates_max = (5.0, 5.0, 5.0)

mesh = TreeMesh(coordinates_min, coordinates_max,initial_refinement_level = 4,n_cells_max = 100_000,periodicity = false)

# Initialize density and pressure perturbation with a Gaussian bump
# that splits into radial waves which are advected with v - c and v + c.
function initial_condition_gauss_wall(x, t, equations::LinearizedEulerEquations3D)
    v1_prime = 0.0
    v2_prime = 0.0
    v3_prime = 0.0
    rho_prime = p_prime = exp(-log(2)/log(ℯ)/9.0*(x[1]^2+x[2]^2+x[3]^2))
    return SVector(rho_prime, v1_prime, v2_prime, v3_prime, p_prime)
end
initial_condition = initial_condition_gauss_wall

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
boundary_conditions = boundary_condition_wall)

###############################################################################
# ODE solvers, callbacks etc.

# At t = 30, the wave moving with v + c crashes into the wall
tspan = (0.0, 24.0)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval = 100)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.9)

save_solution = SaveSolutionCallback(interval = 100, save_initial_solution = true, save_final_solution = true, solution_variables = cons2prim)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution, stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
sol = solve(ode, CarpenterKennedy2N54(thread=OrdinaryDiffEq.True(),williamson_condition = false),
dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
save_everystep = false, callback = callbacks)

# Print the timer summary
summary_callback()

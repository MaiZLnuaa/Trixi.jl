#=
 *  @File          1
 *
 *  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
 *  @Date          Mon Dec 16 2024 3:08:21 PM
 *
 *  @Description   1
 =#

using OrdinaryDiffEq
using Trixi

###############################################################################
# create a restart file

examples_file = "develop_linearizedeuler2d_gauss_tree.jl"
restart_file = "restart_000000288.h5"

trixi_include(@__MODULE__, joinpath(@__DIR__, examples_file))

###############################################################################
# adapt the parameters that have changed compared to "develop_linearizedeuler2d_gauss_tree.jl"

# Note: If you get a restart file from somewhere else, you need to provide
# appropriate setups in the elixir loading a restart file

restart_filename = joinpath("out", restart_file)
mesh = load_mesh(restart_filename)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

tspan = (load_time(restart_filename), 120.0)
dt = load_dt(restart_filename)
ode = semidiscretize(semi, tspan, restart_filename)

# Do not overwrite the initial snapshot written by elixir_advection_extended.jl.
save_solution.condition.save_initial_solution = false

integrator = init(ode, CarpenterKennedy2N54(thread = OrdinaryDiffEq.True(), williamson_condition = false),
                  dt = dt, # solve needs some value here but it will be overwritten by the stepsize_callback
                  save_everystep = false, callback = callbacks);

# Get the last time index and work with that.
load_timestep!(integrator, restart_filename)

###############################################################################
# run the simulation

sol = solve!(integrator)
summary_callback() # print the timer summary
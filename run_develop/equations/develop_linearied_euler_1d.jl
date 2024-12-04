#=
 *  @File          develop_linearied_euler_1d.jl
 *
 *  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
 *  @Date          Tue Dec 03 2024 4:06:04 PM
 *
 *  @Description   write a linearized_euler equation for 1D case
 =#

# must have the following functions:

# @muladd begin
# """
# AbstractLinearizedEulerEquations <: AbstractEquations
# In Trixi, this is a subtype of `AbstractEquations` 
# """
# struct develop_LinearizedEulerEquations1D{RealT <: Real} <:AbstractLinearizedEulerEquations{1, 3}
#     v_mean_global::RealT
#     c_mean_global::RealT
#     rho_mean_global::RealT
# end

# function develop_LinearizedEulerEquations1D(v_mean_global::Real,c_mean_global::Real, rho_mean_global::Real)
#     if rho_mean_global < 0
#         throw(ArgumentError("rho_mean_global must be non-negative"))
#     elseif c_mean_global < 0
#         throw(ArgumentError("c_mean_global must be non-negative"))
#     end
# end

# # Constructor with keywords
# function develop_LinearizedEulerEquations1D(; v_mean_global::Real,c_mean_global::Real,rho_mean_global::Real)
#     return develop_LinearizedEulerEquations1D(v_mean_global, c_mean_global,rho_mean_global)
# end


# """
#     boundary_condition_wall(u_inner, orientation, direction, x, t, surface_flux_function,
#                                 equations::LinearizedEulerEquations1D)

# Boundary conditions for a solid wall.
# """
# function boundary_condition_wall(u_inner, orientation, direction, x, t,
#                                  surface_flux_function,
#                                  equations::develop_LinearizedEulerEquations1D)
#     # Boundary state is equal to the inner state except for the velocity. For boundaries
#     # in the -x/+x direction, we multiply the velocity (in the x direction by) -1.
#     u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3])

#     # Calculate boundary flux
#     if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
#         flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
#     else # u_boundary is "left" of boundary, u_inner is "right" of boundary
#         flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
#     end

#     return flux
# end

# """
# define the physical flux using the calling structure used in Trixi
# """
# @inline function flux()
#     @unpack v_mean_global, c_mean_global, rho_mean_global = equations
#     rho_prime, v1_prime, p_prime = u
#     f1 = v_mean_global * rho_prime + rho_mean_global * v1_prime
#     f2 = v_mean_global * v1_prime + p_prime / rho_mean_global
#     f3 = v_mean_global * p_prime + c_mean_global^2 * rho_mean_global * v1_prime

#     return SVector(f1, f2, f3)
# end

# function varnames(::typeof(cons2cons), ::develop_LinearizedEulerEquations1D)
#     ("rho_prime", "v1_prime", "p_prime")
# end

# # function varnames(::typeof(cons2prim), ::develop_LinearizedEulerEquations1D)
# #     ("rho_prime", "v1_prime", "p_prime")
# # end

# end # @muladd

@muladd begin
    #! format: noindent
    
    @doc raw"""
        develop_LinearizedEulerEquations1D(v_mean_global, c_mean_global, rho_mean_global)
    
    Linearized Euler equations in one space dimension. The equations are given by
    ```math
    \partial_t
    \begin{pmatrix}
        \rho' \\ v_1' \\ p'
    \end{pmatrix}
    +
    \partial_x
    \begin{pmatrix}
        \bar{\rho} v_1' + \bar{v_1} \rho ' \\ \bar{v_1} v_1' + \frac{p'}{\bar{\rho}} \\ \bar{v_1} p' + c^2 \bar{\rho} v_1'
    \end{pmatrix}
    =
    \begin{pmatrix}
        0 \\ 0 \\ 0
    \end{pmatrix}
    ```
    The bar ``\bar{(\cdot)}`` indicates uniform mean flow variables and ``c`` is the speed of sound.
    The unknowns are the perturbation quantities of the acoustic velocity ``v_1'``, the pressure ``p'`` 
    and the density ``\rho'``.
    """
    struct develop_LinearizedEulerEquations1D{RealT <: Real} <:
           AbstractLinearizedEulerEquations{1, 3}
        v_mean_global::RealT
        c_mean_global::RealT
        rho_mean_global::RealT
    end
    
    function develop_LinearizedEulerEquations1D(v_mean_global::Real,
                                        c_mean_global::Real, rho_mean_global::Real)
        if rho_mean_global < 0
            throw(ArgumentError("rho_mean_global must be non-negative"))
        elseif c_mean_global < 0
            throw(ArgumentError("c_mean_global must be non-negative"))
        end
    
        return develop_LinearizedEulerEquations1D(v_mean_global, c_mean_global,
                                          rho_mean_global)
    end
    
    # Constructor with keywords
    function develop_LinearizedEulerEquations1D(; v_mean_global::Real,
                                        c_mean_global::Real, rho_mean_global::Real)
        return develop_LinearizedEulerEquations1D(v_mean_global, c_mean_global,
                                          rho_mean_global)
    end
    
    function varnames(::typeof(cons2cons), ::develop_LinearizedEulerEquations1D)
        ("rho_prime", "v1_prime", "p_prime")
    end
    function varnames(::typeof(cons2prim), ::develop_LinearizedEulerEquations1D)
        ("rho_prime", "v1_prime", "p_prime")
    end
    
    """
        initial_condition_convergence_test(x, t, equations::develop_LinearizedEulerEquations1D)
    
    A smooth initial condition used for convergence tests.
    """
    function initial_condition_convergence_test(x, t, equations::develop_LinearizedEulerEquations1D)
        rho_prime = -cospi(2 * t) * sinpi(2 * x[1])
        v1_prime = sinpi(2 * t) * cospi(2 * x[1])
        p_prime = rho_prime
    
        return SVector(rho_prime, v1_prime, p_prime)
    end
    
    """
        boundary_condition_wall(u_inner, orientation, direction, x, t, surface_flux_function,
                                    equations::develop_LinearizedEulerEquations1D)
    
    Boundary conditions for a solid wall.
    """
    function boundary_condition_wall(u_inner, orientation, direction, x, t,
                                     surface_flux_function,
                                     equations::develop_LinearizedEulerEquations1D)
        # Boundary state is equal to the inner state except for the velocity. For boundaries
        # in the -x/+x direction, we multiply the velocity (in the x direction by) -1.
        u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3])
    
        # Calculate boundary flux
        if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
            flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
        else # u_boundary is "left" of boundary, u_inner is "right" of boundary
            flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
        end
    
        return flux
    end
    
    # Calculate 1D flux for a single point
    @inline function flux(u, orientation::Integer, equations::develop_LinearizedEulerEquations1D)
        @unpack v_mean_global, c_mean_global, rho_mean_global = equations
        rho_prime, v1_prime, p_prime = u
        f1 = v_mean_global * rho_prime + rho_mean_global * v1_prime
        f2 = v_mean_global * v1_prime + p_prime / rho_mean_global
        f3 = v_mean_global * p_prime + c_mean_global^2 * rho_mean_global * v1_prime
    
        return SVector(f1, f2, f3)
    end
    
    @inline have_constant_speed(::develop_LinearizedEulerEquations1D) = True()
    
    @inline function max_abs_speeds(equations::develop_LinearizedEulerEquations1D)
        @unpack v_mean_global, c_mean_global = equations
        return abs(v_mean_global) + c_mean_global
    end
    
    @inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                         equations::develop_LinearizedEulerEquations1D)
        @unpack v_mean_global, c_mean_global = equations
        return abs(v_mean_global) + c_mean_global
    end
    
    # Calculate estimate for minimum and maximum wave speeds for HLL-type fluxes
    @inline function min_max_speed_naive(u_ll, u_rr, orientation::Integer,
                                         equations::develop_LinearizedEulerEquations1D)
        min_max_speed_davis(u_ll, u_rr, orientation, equations)
    end
    
    # More refined estimates for minimum and maximum wave speeds for HLL-type fluxes
    @inline function min_max_speed_davis(u_ll, u_rr, orientation::Integer,
                                         equations::develop_LinearizedEulerEquations1D)
        @unpack v_mean_global, c_mean_global = equations
    
        λ_min = v_mean_global - c_mean_global
        λ_max = v_mean_global + c_mean_global
    
        return λ_min, λ_max
    end
    
    # Convert conservative variables to primitive
    @inline cons2prim(u, equations::develop_LinearizedEulerEquations1D) = u
    @inline cons2entropy(u, ::develop_LinearizedEulerEquations1D) = u
    end # muladd
    
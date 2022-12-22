module DifferentiableOSQP


using LinearAlgebra, OSQP, SparseArrays, ForwardDiff

include("solve_wrappers.jl")
include("derivatives.jl")

export solve, solve_and_jac



end

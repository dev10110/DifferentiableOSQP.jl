# DifferentiableOSQP

<!--- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dev10110.github.io/DifferentiableOSQP.jl/stable/) --->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dev10110.github.io/DifferentiableOSQP.jl/dev/)
[![Build Status](https://github.com/dev10110/DifferentiableOSQP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dev10110/DifferentiableOSQP.jl/actions/workflows/CI.yml?query=branch%3Amain)



This package provides a thin wrapper of [OSQP.jl](https://github.com/osqp/OSQP.jl), but also provides the ability to differentiate through the quadratic progam, based on the equations in [OptNet](https://arxiv.org/abs/1703.00443). 

The package exports 2 commands: `solve` and `solve_and_jac`.


## Installation

Activate your environment and simply add `DifferentialOSQP`:
```
] add https://github.com/dev10110/DifferentiableOSQP.jl
```


## Interface 1: Solve a QP
```
   x = solve(P, q, A, u; kwargs...)
```
Solves the quadratic program:
```math
\begin{aligned}
  \min_x   &\quad \frac{1}{2} x^T P x + q^T x \\ 
  s.t.    &\quad A x \leq u
\end{aligned}
```
where `kwargs` are keyword arguments passed into `OSQP.setup!`

To include equality constraints, for example ``G x = h``, modify `A, u` matrices as:
```
A = [A ; G ; -G]
u = [u ; h; -h]
```
This introduces the constraints ``G x \leq h`` and ``G x \geq h``, allowing equalities to be handled. Yes, this is rather inefficient, but the easiest way to solve the problem I think.


## Interface 2: Solve a Parameteric QP

```
  x = solve(θ, P_fn, q_fn, A_fn, u_fn; kwargs)
```
Solves the quadratic program:
```math
\begin{aligned}
  \min_x   &\quad \frac{1}{2} x^T P(\theta) x + q(\theta)^T x \\ 
  s.t.    &\quad A(\theta) x \leq u(\theta)
\end{aligned}
```
i.e., assumes the `P_fn, q_fn, A_fn, u_fn` are functions of `θ`. Note, the shape and size of each output must be correct - for example, `A_fn(θ)` must return a matrix, and `q_fn(θ)` must return a vector. 

## Interface 3: Jacobians of a QP

Thinking about interface 2, notice that a QP solver is essentially a function
```math
QP : \mathbb{R}^p \to \mathbb{R}^n\\
x = QP(θ)
```

Therefore, the Jacobian of the QP 
```math
J = \frac{\partial x}{\partial \theta}
```

If we want the jacobian of ``QP``, we can call it as the following:
```
x, J = solve_and_jac(θ, P_fn, q_fn, A_fn, u_fn; kwargs)
```
which gives the optimal solution `x`, and the jacobian `J`


## Interaface 4: `ForwardDiff`

For convenience, we overloaded `solve` to handle Dual numbers. This means we can directly use `ForwardDiff.jacobian`, as in the following example:

```
J = ForwardDiff.jacobian(θ -> solve(θ, P_fn, q_fn, A_fn, u_fn), θ0)
```

or more explictly
```
function parameteric_qp(θ)
  return solve(θ, P_fn, q_fn, A_fn, u_fn)
end

J = ForwardDiff.jacobian(parametric_qp, θ0)
```


## Warning

Naturally, not all QPs are differentiable. This library will always return a derivative, but doesnt check/warn if the derivative doesnt exist. I want to add this functionality in the future. See [this](https://doi.org/10.1109/CDC.2013.6760327) or [this](https://doi.org/10.1016/0022-247X(67)90163-1) paper for some results on existence of derivatives/Lipschitz continuity of QPs. 

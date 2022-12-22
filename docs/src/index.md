```@meta
CurrentModule = DifferentiableOSQP
```

# DifferentiableOSQP

Documentation for [DifferentiableOSQP](https://github.com/dev10110/DifferentiableOSQP.jl).

```@index
```

This package provides a thin wrapper of [OSQP.jl](https://github.com/osqp/OSQP.jl), but also provides the ability to differentiate through the quadratic progam, based on the equations in [OptNet](https://arxiv.org/abs/1703.00443). 

The package exports 2 commands: `solve` and `solve_and_jac`.


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




```@autodocs
Modules = [DifferentiableOSQP]
```


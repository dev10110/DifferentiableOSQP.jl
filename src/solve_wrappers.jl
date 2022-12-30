

function solve(P, q, A, u; verbose=false, eps_rel=1e-05, eps_abs=1e-5, kwargs...)
  """
  Uses OSQP to solve the QP:
  
  min_x 1/2 x' P x + q' x
  s.t.  A ≤ u

  To include equality constraints, extend A such that A <= u and -A <= -u are constraints.
  
  """
 
  m = length(u)
  l = fill(-Inf, m)

  # setup OSQP problem
  prob = OSQP.Model()
  OSQP.setup!(prob; P=sparse(P), q=q, A=sparse(A), l=l, u=u, 
              verbose=verbose,
              eps_rel=eps_rel,
              eps_abs=eps_abs,
              kwargs...)
  results = OSQP.solve!(prob)
  
  # check the result
  @assert results.info.status in OSQP.SOLUTION_PRESENT "OSQP could not find a solution to the QP"
  if results.info.status != :Solved
      @warn "Solution is not optimal. Status: $(results.info.status)"
  end
  
  # return the primal
  return results.x
  
end

function solve(θ, P_fn, q_fn, A_fn, u_fn; verbose=false,  kwargs...)
  """
  Uses OSQP to solve the QP:
  
  min_x 1/2 x' P(θ) x + q(θ)' x
  s.t.  A(θ) ≤ u(θ)

  To include equality constraints, extend A such that A <= u and -A <= -u are constraints.
  
  """
  
  return solve(P_fn(θ), q_fn(θ), A_fn(θ), u_fn(θ); verbose=verbose, kwargs...)
  
end


function solve(θ::Vector{ForwardDiff.Dual{T, V, N}}, P_fn, q_fn, A_fn, u_fn; verbose=false,  kwargs...)  where {T, V, N}
  "specialization to handle forward diff requests"
  
  p = length(θ)
  
  # get just the values
  θ_vals = [θi.value for θi in θ]
  
  # get the jacobian of θ
  Jθ = zeros(V, p, N)
  for i=1:p, j=1:N
      Jθ[i, j] = θ[i].partials[j]
  end
  
  # actually solve the QP, and get its derivative
  x, J_QP = solve_and_jac(θ_vals, P_fn, q_fn, A_fn, u_fn; verbose=verbose, kwargs...)

  # construct the overall jacobian
  J = J_QP * Jθ
  
  # construct the dual output
  dual_output = [ForwardDiff.Dual{T}(x[i], J[i, :]...) for i in eachindex(x)]
  
  return dual_output
  
end

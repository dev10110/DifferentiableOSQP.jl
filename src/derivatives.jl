

function primal_dual(P, q, A, u; verbose=false, kwargs...)
  """
  Uses OSQP to solve the QP:
  
  min_x 1/2 x' P x + q' x
  s.t.  A ≤ u

  Returns (primal, dual) as a tuple
  
  """
 
  m = length(u)
  l = fill(-Inf, m)

  # setup OSQP problem
  prob = OSQP.Model()
  OSQP.setup!(prob; P=sparse(P), q=q, A=sparse(A), l=l, u=u, verbose=verbose, kwargs...)
  results = OSQP.solve!(prob)
  
  # check the result
  @assert results.info.status in OSQP.SOLUTION_PRESENT "OSQP could not find a solution to the QP"
  if results.info.status != :Solved
      @warn "Solution is not optimal. Status: $(results.info.status)"
  end
  
  # return the primal
  return results.x, results.y
  
end


# implements equation 6 of the paper, but without equality constraints
function qp_deriv(x_s, y_s, P, q, A, u, dP, dq, dA, du)
  
  n = length(x_s)
  
  Dy_s = Diagonal(y_s)
  
  M = [[P ;; A']; [Dy_s * A ;; Diagonal(A*x_s - u)]]
  
  N = [
       dP * x_s + dq + dA' * y_s;
       Dy_s * dA * x_s - Dy_s*du;
      ]
  
  # let dl = [dx, dy, dν]
  dl = - M \ N
  
  dx = dl[1:n]
  #dy = dl[(n+1):(n+p)]
        
  return dx
  
end


function solve_and_jac(θ, P_fn, q_fn, A_fn, u_fn; verbose=false, kwargs...)
        
  P = P_fn(θ)
  q = q_fn(θ)
  A = A_fn(θ)
  u = u_fn(θ)
  
  m, n = size(A)
  p = length(θ)
  
  res_x, res_y = primal_dual(P, q, A, u; verbose=false, kwargs...)
  
  # get jacobians
  JP = ForwardDiff.jacobian(P_fn, θ)
  JA = ForwardDiff.jacobian(A_fn, θ)
  Jq = ForwardDiff.jacobian(q_fn, θ)
  Ju = ForwardDiff.jacobian(u_fn, θ)
  
  # allocate memory for the J_QP
  J_QP = zeros(n, p)

  # for each of the dimensions of θ
  for i=1:p
      dP = reshape(JP[:,i], n,n)
      dA = reshape(JA[:,i], m,n)
      dq = Jq[:,i]
      du = Ju[:,i]
      
      J_QP[:, i] = qp_deriv(res_x, res_y, P, q, A, u, dP, dq, dA, du)
  end
  
  return res_x, J_QP
  
end

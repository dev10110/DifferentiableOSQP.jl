using LinearAlgebra, Test, DifferentiableOSQP, ForwardDiff, FiniteDifferences

@testset "solve a qp" begin

  P = 1.0*I(2)

  A = [[1;; 0]; [1;; 1.]]

  q = [0,0.]

  u = [-1, -2.]


  x = solve(P, q, A, u)

  @test norm(x  - [-1, -1.]) <= 1e-3


end

@testset "parameteric qp" begin

  function gen_P(θ)
    P = diagm([θ[1], θ[2]])
  end

  function gen_q(θ)
      q = [0,0.]
  end

  function gen_A(θ)
      A = [[θ[3];; 0];]
  end

  function gen_u(θ)    
      u = [θ[4]]
  end

  th = [1,1,1,1.]

  x = solve(th, gen_P, gen_q, gen_A, gen_u)

  @test norm(x  - [0, 0.]) <= 1e-4


  th = [1,1,1,-1.]

  x = solve(th, gen_P, gen_q, gen_A, gen_u)

  @test norm(x  - [-1, 0]) <= 1e-3


end

@testset "differentiate qp - A, u" begin

  function gen_P(θ)
    P = diagm([θ[1], θ[2]])
  end

  function gen_q(θ)
      q = [0,0.]
  end

  function gen_A(θ)
      A = [[θ[3];; 0];]
  end

  function gen_u(θ)    
      u = [θ[4]]
  end

  
  function parametric_qp(th)
    return solve(th, gen_P, gen_q, gen_A, gen_u)
  end

  th = [1,1,1,-1.]
  J = ForwardDiff.jacobian(parametric_qp, th)
  J_analytic = [[0 ;; 0;; -th[4]/th[3]^2 ;; 1.0/th[3]] ; [0 ;; 0 ;; 0;; 0]]

  @test norm(J - J_analytic) <= 1e-3


  th = [1,2,3,-4.]
  J = ForwardDiff.jacobian(parametric_qp, th)
  J_analytic = [[0 ;; 0;; -th[4]/th[3]^2 ;; 1.0/th[3]] ; [0 ;; 0 ;; 0;; 0]]

  @test norm(J - J_analytic) <= 1e-3




end


@testset "differentiate qp - q" begin

  function gen_P(θ)
    P = diagm([1, 1.])
  end

  function gen_q(θ)
      q = [θ[1],θ[2]]
  end

  function gen_A(θ)
      A = [[1;; 1.];]
  end

  function gen_u(θ)    
      u = [0.0]
  end

  
  function parametric_qp(th)
    return solve(th, gen_P, gen_q, gen_A, gen_u)
  end

  th = [1,1.]
  J = ForwardDiff.jacobian(parametric_qp, th)
  J_analytic = [[-1 ;; 0] ; [0 ;; -1]]

  @test norm(J - J_analytic) <= 1e-3


  th = [-1,-2.]
  J = ForwardDiff.jacobian(parametric_qp, th)
  J_analytic = [[-0.5 ;; 0.5] ; [0.5 ;; -0.5]]

  @test norm(J - J_analytic) <= 1e-3


end

@testset "differentiate qp - P" begin

  function gen_P(θ)
    P = [[θ[1] ;; -2.0] ; [-2 ;; θ[2]]]
  end

  function gen_q(θ)
      q = [1, 1.]
  end

  function gen_A(θ)
      A = [[1;; +1.];]
  end

  function gen_u(θ)
      u = [0.0]
  end

  
  function parametric_qp(th)
    return solve(th, gen_P, gen_q, gen_A, gen_u)
  end

  # analytic solutions derived using mathematica
  # for solutions to exist, need th2 > 0 and th1 th2 > 4
  q1 = 5.0
  q2 = 1.0
  th = [q1, q2]
  J = ForwardDiff.jacobian(parametric_qp, th)

  
  J_analytic = [ [(q2*(2 + q2))/(-4 + q1*q2)^2 ;;  (2*(2 + q1))/(-4 + q1*q2)^2] ;
                 [(2*(2 + q2))/(-4 + q1*q2)^2 ;; (q1*(2 + q1))/(-4 + q1*q2)^2 ]
                 ]


  @test norm(J - J_analytic) <= 1e-4


   # for solutions to exist, need th2 > 0 and th1 th2 > 4
   q1 = 24.0
   q2 = 0.5
   th = [q1, q2]
   x, J = solve_and_jac(th, gen_P, gen_q, gen_A, gen_u)

   x_analytic = [(2 + q2)/(4 - q1*q2), (2 + q1)/(4 - q1*q2)]
   
   J_analytic = [ [(q2*(2 + q2))/(-4 + q1*q2)^2 ;;  (2*(2 + q1))/(-4 + q1*q2)^2] ;
                  [(2*(2 + q2))/(-4 + q1*q2)^2 ;; (q1*(2 + q1))/(-4 + q1*q2)^2 ]
                  ]
 
 
  @test norm(x - x_analytic) <= 1e-3
   @test norm(J - J_analytic) <= 1e-3

end



# @testset "differentiate qp - finitediff" begin

#   function gen_P(θ)
#     M = diagm([θ[1], θ[2], θ[2]*θ[1]]) 
#     P = M * M' + 1e-2*I
#   end

#   function gen_q(θ)
#       q = [θ[2],θ[2]^2, θ[1]]
#   end

#   function gen_A(θ)
#       A = [[θ[1];; sin(θ[2]) + exp(θ[1]) ;; 3.0];]
#   end

#   function gen_u(θ)    
#       u = [θ[2]]
#   end

  
#   function parametric_qp(th)
#     return solve(th, gen_P, gen_q, gen_A, gen_u)
#   end

#   for i=1:15
#     th = randn(2)
#     J = ForwardDiff.jacobian(parametric_qp, th)
#     J_FD= FiniteDifferences.jacobian(central_fdm(12, 1), parametric_qp, th)[1]

#     @test norm(J - J_FD) <= 1e-2

#   end





# end

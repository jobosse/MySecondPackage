using DifferentialEquations, Random, Test

include("../src/Predefined_Zero_D_EBM.jl")

function zero_d_ebm(u,p,t) 
    α, S₀, σ = p 
    return (1 - α)*(S₀/4) - σ*u[1]^4
end
p = [rand(), rand()*2000, rand()*1e-8]
u0 = rand()*200
t = rand()*2
prob = ODEProblem(zero_d_ebm, u0, t, p)
sol_1 = solve(prob)
sol_2 = solve(SimplestModel(;u0=u0,p=p,T=t))
@test all(sol_1 .≈ sol_2)

function fixed_point(p,ε)
    α, S₀, σ = p
    return (((1-α)*S₀)/(4*ε*σ))^(1/4)
end

sol_3 = solve(SimplestModel(;u0=u0, p=p, T=1000.))
@test abs(sol_3.u[end] - fixed_point(p,1)) < 0.1
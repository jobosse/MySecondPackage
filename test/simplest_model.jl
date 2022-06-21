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
sol = solve(prob)
@test all(sol .≈ solve(SimplestModel(;u0=u0,p=p,T=t)))
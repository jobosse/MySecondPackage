using DifferentialEquations
"""
SimplestModel(u0=273.,p=[0.3, 1363, 5.67e-8], T=2.5)

Returns DifferentialEquations.ODEProblem of the simplest zero dimensional energy model.

---
Arguments:
- u0::Float64, initial temperature
- p::Array{Float64}: [α, S₀, σ] = p
- t:: Float64, end time of the integration

---
Returns:
- DifferentialEquations.ODEProblem object.
"""
function SimplestModel(;u0=273.,p=[0.3, 1363, 5.67e-8], T=2.5)
    function zero_d_ebm(u,p,t) 
        α, S₀, σ = p 
        return (1 - α)*(S₀/4) - σ*u[1]^4
    end 
    prob = ODEProblem(zero_d_ebm, u0, (0.,T),p)
    return prob
end




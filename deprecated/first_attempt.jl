################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# CoolingInfrastructureModels.jl                                               #
# A Julia package for cooling infrastructure optimization within a large       #
# industrial site.                                                             #
# See http://github.com/timmyfaraday/CoolingInfrastructureOptimization.jl      #
################################################################################

# In this file, the feasible space is modeled for a cooling infrastructure using
# brackish water. It consists of two water networks: 
# - a 'cold' pressurized network with one source
# - a 'hot' gravitational network with one drain
# Furthermore, it has one fixed unit, i.e., a plant, and one dispatchable unit:
# a cooling tower. 

# Notes: positive pipe flow is defined from the lowest node number to the highest
#        positive pump flow is defined based on its direction

using JuMP

# sets
## nodes
𝓝       = 1:11

## pipes
𝓐       = Set(["a₁","a₂","a₃","a₄","a₅","a₆","a₇","a₈","a₉","a₁₀","a₁₁"])
𝓐ᵗᵖˡ    = Set([("a₁",2,3),("a₂",2,4),("a₃",2,5),("a₄",3,5),("a₅",4,5),("a₆",6,7),
                ("a₇",6,8),("a₈",6,9),("a₉",7,9),("a₁₀",8,9),("a₁₁",8,10)])

## pumps
𝓟       = Set(["p₁","p₂"])
𝓟ᵗᵖˡ    = Set([("p₁",1,3),("p₂",11,4)])

## edges
𝓔       = union(𝓐,𝓟)

## dispatchable units
𝓤ᵈ      = Set(["c₁"])
𝓤ᵈ⁻ᵗᵖˡ  = Set([("c₁",6,11)])

## fixed units
𝓤ᶠ      = Set(["f₁"])
𝓤ᶠ⁻ᵗᵖˡ  = Set([("f₁",5,7)])

## units
𝓤       = union(𝓤ᵈ,𝓤ᶠ)

## sources
𝓢       = Set(["r₁","r₂"])
𝓢ᵗᵖˡ    = Set([("r₁",1),("r₂",10)])

# parameters
## head
hˢʳᶜ    = Dict( "r₁" => 10.0,
                "r₂" => 5.0)
hᵐⁱⁿ     = Dict( 1 => 0.0)
hᵐᵃˣ    = Dict( 1 => 300.0)
## volumetric flow
qᶠ      = Dict( "f₁" => 100.0)
qᵐⁱⁿ     = Dict( "c₁" => 0.0,
                "p₁" => 0.0)
qᵐᵃˣ    = Dict( "c₁" => 200.0,
                "p₁" => 200.0)
vᵐⁱⁿ     = Dict( "a₁" => 0.5)
vᵐᵃˣ    = Dict( "a₁" => 3.0)
vⁱⁿᶠ     = Dict( "a₁" => 5.0)

## pipe parameters
d       = Dict( "a₁" => 2.0)
l       = Dict( "a₁" => 180.0)
rˡ      = Dict( "a₁" => 3.0144e-5)
rᶜ      = Dict( "a₁" => 0.1703)

## pump parameters
Δh      = Dict( "p₁" => 200.0)
cᵃ      = Dict( "p₁" => 0.010)

# model
model = Model()

# variables
@variable(model, h[i ∈ 𝓝])
@variable(model, 0.0 <= q⁻[a ∈ 𝓐])
@variable(model, 0.0 <= q⁺[a ∈ 𝓐])
@variable(model, q[e ∈ 𝓔])
@variable(model, q[s ∈ 𝓢])
@variable(model, q[u ∈ 𝓤])

# constraints
## head constraints
@constraint(model,  [(s,i) ∈ 𝓢ᵗᵖˡ, k ∈ 𝓚], 
                    h[i,k] == hˢʳᶜ[s])
@constraint(model,  [i ∈ 𝓝], 
                    [(hᵐᵃˣ[i] - mean(h[i,:],mop)) / λʰ[i]; variancex(h[i,:],mop)] in SecondOrderCone())
## volumetric flow rate
@constraint(model,  [u ∈ 𝓤ᵈ, k ∈ 𝓚], 
                    qᵐⁱⁿ[u] <= q[u,k])
@constraint(model,  [u ∈ 𝓤ᵈ, k ∈ 𝓚], 
                    q[u,k] <= qᵐᵃˣ[u])
@constraint(model,  [u ∈ 𝓤ᶠ, k ∈ 𝓚], 
                    q[u,k] == qᶠ[u])
@constraint(model,  [a ∈ 𝓐, k ∈ 𝓚], 
                    q[a,k] == q⁺[a,k] - q⁻[a,k])
@constraint(model,  [a ∈ 𝓐, k ∈ 𝓚], 
                    q⁻[a,k] <= (1 - x[a,k]) * π / 4 * d[a]^2 * vⁱⁿᶠ[a])
@constraint(model,  [a ∈ 𝓐, k ∈ 𝓚], 
                    q⁺[a,k] <= x[a,k] * π / 4 * d[a]^2 * vⁱⁿᶠ[a])
@constraint(model,  [i ∈ 𝓝], 
                    [(mean(h[i,:],mop) - hᵐⁱⁿ[i]) / λʰ[i]; variancex(h[i,:],mop)] in SecondOrderCone())
@constraint(mode,   [p ∈ 𝓟],
                    [(qᵐᵃˣ[p] - mean(q[p,:],mop)) / λᵖ[p]; variancex(q[p,:],mop)] in SecondOrderCone())
@constraint(mode,   [p ∈ 𝓟],
                    [(mean(q[p,:],mop) - qᵐⁱⁿ[p]) / λᵖ[p]; variancex(q[p,:],mop)] in SecondOrderCone())
@constraint(model,  [a ∈ 𝓐],
                    [(π / 4 * d[a]^2 * vᵐᵃˣ[a] - mean(q⁻[a,:],mop)) / λᵛ[a]; variancex(q⁻[a,:],mop)] in SecondOrderCone())
@constraint(model,  [a ∈ 𝓐],
                    [(mean(q⁻[a,:],mop) - π / 4 * d[a]^2 * vᵐⁱⁿ[a]) / λᵛ[a]; variancex(q⁻[a,:],mop)] in SecondOrderCone())
@constraint(model,  [a ∈ 𝓐],
                    [(π / 4 * d[a]^2 * vᵐᵃˣ[a] - mean(q⁺[a,:],mop)) / λᵛ[a]; variancex(q⁺[a,:],mop)] in SecondOrderCone())
@constraint(model,  [a ∈ 𝓐],
                    [(mean(q⁺[a,:],mop) - π / 4 * d[a]^2 * vᵐⁱⁿ[a]) / λᵛ[a]; variancex(q⁺[a,:],mop)] in SecondOrderCone())
## potential-flow coupling constraints
@constraint(model,  [(p,i,j) ∈ 𝓟ᵗᵖˡ, k ∈ 𝓚], 
                    h[j,k] - h[i,k] == Δh[p] - cᵃ[p] * q[p,k])
@constraint(model,  [(a,i,j) ∈ 𝓐ᵗᵖˡ, k ∈ 𝓚], 
                    h[i,k] - h[j,k] == rˡ[a] * l[a] * q[a,k] + rᶜ[a] * (1 - 2 * x[a,k]))
## mass conservation constraint
@constraint(model,  [i ∈ 𝓝, k ∈ 𝓚], 
                    sum(q[e,k] for e ∈ 𝓔⁻[i,k]) - sum(q[e,k] for e ∈ 𝓔⁺[i,k]) +
                    sum(q[u,k] for u ∈ 𝓤⁻[i,k]) - sum(q[u,k] for u ∈ 𝓤⁺[i,k]) +
                    sum(q[s,k] for s ∈ 𝓢ⁱ[i,k]))
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

##### NEED TO BE ADJUSTED FOR THE DETERMINISTIC LOOSE PARAMETERS!!!!!

# using pkgs
using Ipopt, JuMP, PolyChaos

# sets
N = 1:11
E = Set(["a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","p1","p2"])
S = Set(["s1","s2"])
U = Set(["u1","u2"])

# subset
A    = Set(["a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11"])
Aᵗᵖˡ = Set([("a1",2,3),("a2",2,4),("a3",2,5),("a4",3,5),("a5",4,5),("a6",6,7),
            ("a7",6,8),("a8",6,9),("a9",7,9),("a10",8,9),("a11",8,10)])
E⁻    = Dict(1 => ["p1"],
            2 => ["a1","a2","a3"],
            3 => ["a4"],
            4 => ["a5"],
            5 => [],
            6 => ["a6","a7","a8"],
            7 => ["a9"],
            8 => ["a10","a11"],
            9 => [],
            10 => [],
            11 => ["p2"])
E⁺    = Dict(1 => [],
            2 => [],
            3 => ["p1","a1"],
            4 => ["p2","a2"],
            5 => ["a3","a4","a5"],
            6 => [],
            7 => ["a6"],
            8 => ["a7"],
            9 => ["a8","a9","a10"],
            10 => ["a11"],
            11 => [])
P    = Set(["p1","p2"])
Pᵗᵖˡ = Set([("p1",1,3),("p2",11,4)])
Sⁱ    = Dict(1 => ["s1"],
            2 => [],
            3 => [],
            4 => [],
            5 => [],
            6 => [],
            7 => [],
            8 => [],
            9 => [],
            10 => ["s2"],
            11 => [])
Sᵗᵖˡ = Set([("s1",1),("s2",10)])
Uᵈ   = Set(["u1"])
Uᶠ   = Set(["u2"])
U⁻   = Dict(1 => [],
            2 => [],
            3 => [],
            4 => [],
            5 => ["u2"],
            6 => ["u1"],
            7 => [],
            8 => [],
            9 => [],
            10 => [],
            11 => [])
U⁺   = Dict(1 => [],
            2 => [],
            3 => [],
            4 => [],
            5 => [],
            6 => [],
            7 => ["u2"],
            8 => [],
            9 => [],
            10 => [],
            11 => ["u1"]) 
X    = union(E,S,U)

# stochastic sets
deg = 1
opq = [GaussOrthoPoly(deg; Nrec=5*deg)]
mop = MultiOrthoPoly(opq, deg)
K   = 1:mop.dim
λʰ⁻  = 1.6*ones(length(N))
λʰ⁺  = 1.6*ones(length(N))
λᵖ⁻  = Dict(p => 1.6 for p in P)
λᵖ⁺  = Dict(p => 1.6 for p in P)
λᵛ⁻  = Dict(a => 1.6 for a in A)
λᵛ⁺  = Dict(a => 1.6 for a in A)
T2  = Tensor(2,mop)
T3  = Tensor(3,mop)

# parameters
## head parameters
hˢʳᶜ = Dict("s1" => 4.3,
            "s2" => 2.0)
hᵐⁱⁿ  = Dict(i => 0 for i in N)
hᵐᵃˣ = Dict(i => 1500 for i in N)
## (volumetric) flow rate parameters
qᶠ   = Dict("u2" => convert2affinePCE(10.0,1.0,opq[1]))
qᵐⁱⁿ  = Dict(x => 0.0 for x in union(P,Uᵈ))
qᵐᵃˣ = Dict(x => 15.0 for x in union(P,Uᵈ))
vᵐⁱⁿ = Dict(a => 0.0 for a in A)
vᵐᵃˣ = Dict(a => 5.0 for a in A)
## pump parameters
Δh    = Dict(p => 40.0 for p in P)
cᵃ    = Dict(p => 1.75 for p in P)
## pipe parameters
d     = Dict(a => 2.0 for a in A)
r     = Dict(a => 6.3967e-5 for a in A)
l     = Dict(a => 100.0 for a in A)

# functions
function _f()
    return function(x::Float64)
        return sign(x)
    end
end
function _df()
    return function(x::Float64)
        return sign(x)
    end
end
function _d2f()
    return function(x::Float64)
        return x != 0.0 ? sign(x) : prevfloat(Inf)
    end
end
function signed_volumetric_flow_rate_args()
    return (:signed_volumetric_flow_rate, 1, _f(), _df(), _d2f())
end

# model
model = Model(Ipopt.Optimizer)

# register signed volumetric flow rate
register(model, signed_volumetric_flow_rate_args()...)

# variables
@variable(  model,
            h[i in N, k in K])
@variable(  model,
            q[x in X, k in K])
@variable(  model,
            q²[e in E, k in K])

# constraints
## head constraints
@constraint(model,
            [(s,i) in Sᵗᵖˡ, k in K],
            h[i,k] == hˢʳᶜ[s])
@constraint(model,
            [i in N],
            var(h[i,:],T2) <= ((mean(h[i,:],mop) - hᵐⁱⁿ[i]) / λʰ⁻[i])^2)
@constraint(model,
            [i in N],
            var(h[i,:],T2) <= ((hᵐᵃˣ[i] - mean(h[i,:],mop)) / λʰ⁺[i])^2)
## volumetric flow rate constraints
@constraint(model,
            [u in Uᵈ, k in K],
            qᵐⁱⁿ[u] <= q[u,k])
@constraint(model,
            [u in Uᵈ, k in K],
            q[u,k] <= qᵐᵃˣ[u])
@constraint(model,
            [u in Uᶠ, k in K],
            q[u,k] == qᶠ[u][k])
@constraint(model,
            [e in E, k in K],
            q²[e,k] == sum(T3.get([k-1,k₁-1,k₂-1]) * q[e,k₁] * q[e,k₂] for k₁ in K, k₂ in K))
@constraint(model,
            [p in P],
            var(q²[p,:],T2) <= ((mean(q²[p,:],mop) - (qᵐⁱⁿ[p])^2) / λᵖ⁻[p])^2)
@constraint(model,
            [p in P],
            var(q²[p,:],T2) <= (((qᵐᵃˣ[p])^2 - mean(q²[p,:],mop)) / λᵖ⁺[p])^2)
@constraint(model,
            [a in A],
            var(q²[a,:],T2) <= ((mean(q²[a,:],mop) - (π / 4 * d[a]^2 * vᵐⁱⁿ[a])^2) / λᵛ⁻[a])^2)
@constraint(model,
            [a in A],
            var(q²[a,:],T2) <= (((π / 4 * d[a]^2 * vᵐᵃˣ[a])^2 - mean(q²[a,:],mop)) / λᵛ⁺[a])^2)
## potential flow coupling constraints
@constraint(model,
            [(p,i,j) in Pᵗᵖˡ, k in K],
            h[j,k] - h[i,k] == Δh[p] - cᵃ[p] * q²[p,k])
@NLconstraint(model,
            [(a,i,j) in Aᵗᵖˡ, k in K],
            h[i,k] - h[j,k] == r[a] * l[a] * q²[a,k] * signed_volumetric_flow_rate(q[a,k]))
## mass conservation constraint
@constraint(model,
            [i in N, k in K],
            sum(q[e,k] for e in E⁺[i]) - sum(q[e,k] for e in E⁻[i]) +
            sum(q[u,k] for u in U⁺[i]) - sum(q[u,k] for u in U⁻[i]) + 
            sum(q[s,k] for s in Sⁱ[i]) == 0)

# objective
# @objective(model, Max, q["u1",1])
# @objective(model, Max, q["u1"]) --> feasibility problem, no issue, giving it a objective seems to cause problems, needs further investigation.

# optimize
optimize!(model)
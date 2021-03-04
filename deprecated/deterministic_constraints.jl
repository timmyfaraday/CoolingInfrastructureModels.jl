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

###### FULLY FIGURED OUT - I THINK

# using pkgs
using Ipopt, JuMP

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

# parameters
## head parameters
hˢʳᶜ = Dict("s1" => 5.0,
            "s2" => 2.0)
hᵐⁱⁿ  = Dict(i => 0 for i in N)
hᵐᵃˣ = Dict(i => 1000 for i in N)
## (volumetric) flow rate parameters
qᶠ   = Dict("u2" => 10.0)
qᵐⁱⁿ  = Dict(x => 0.0 for x in union(P,Uᵈ))
qᵐᵃˣ = Dict(x => 12.0 for x in union(P,Uᵈ))
vᵐⁱⁿ = Dict(a => 0.0 for a in A)
vᵐᵃˣ = Dict(a => 5.0 for a in A)
## pump parameters
Δh    = Dict(p => 1000.0 for p in P)
cᵃ    = Dict(p => 0.10 for p in P)
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
            h[i in N])
@variable(  model,
            q[x in X])
@variable(  model,
            q²[e in E])

# constraints
## head constraints
@constraint(model,
            [(s,i) in Sᵗᵖˡ],
            h[i] == hˢʳᶜ[s])
@constraint(model,
            [i in N],
            hᵐⁱⁿ[i] <= h[i])
@constraint(model,
            [i in N],
            h[i] <= hᵐᵃˣ[i])
## volumetric flow rate constraints
@constraint(model,
            [u in Uᵈ],
            qᵐⁱⁿ[u] <= q[u])
@constraint(model,
            [u in Uᵈ],
            q[u] <= qᵐᵃˣ[u])
@constraint(model,
            [u in Uᶠ],
            q[u] == qᶠ[u])
@constraint(model,
            [e in E],
            q²[e] == q[e] * q[e])
@constraint(model,
            [p in P],
            (qᵐⁱⁿ[p])^2 <= q²[p])
@constraint(model,
            [p in P],
            q²[p] <= (qᵐᵃˣ[p])^2)
@constraint(model,
            [a in A],
            (π / 4 * d[a]^2 * vᵐⁱⁿ[a])^2 <= q²[a])
@constraint(model,
            [a in A],
            q²[a] <= (π / 4 * d[a]^2 * vᵐᵃˣ[a])^2)
## potential flow coupling constraints
@constraint(model,
            [(p,i,j) in Pᵗᵖˡ],
            h[j] - h[i] == Δh[p] - cᵃ[p] * q²[p])
@NLconstraint(model,
            [(a,i,j) in Aᵗᵖˡ],
            h[i] - h[j] == r[a] * l[a] * q²[a] * signed_volumetric_flow_rate(q[a]))
## mass conservation constraint
@constraint(model,
            [i in N],
            sum(q[e] for e in E⁺[i]) - sum(q[e] for e in E⁻[i]) +
            sum(q[u] for u in U⁺[i]) - sum(q[u] for u in U⁻[i]) + 
            sum(q[s] for s in Sⁱ[i]) == 0)

# objective
# @objective(model, Max, q["u1"]) --> feasibility problem, no issue, giving it a objective seems to cause problems, needs further investigation.

# optimize
optimize!(model)
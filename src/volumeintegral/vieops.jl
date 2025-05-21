abstract type VIEOperator <: IntegralOperator end

abstract type VolumeOperator <: VIEOperator end #6D integral: ∫∫∫_Ω  ∫∫∫_Ω

abstract type BoundaryOperator <: VIEOperator end
abstract type BoundaryOperator_ΓΩ <: BoundaryOperator end #5D integral: ∫∫_Γ ∫∫∫_Ω
abstract type BoundaryOperator_ΩΓ <: BoundaryOperator end #5D integral: ∫∫∫_Ω ∫∫_Γ 

struct VIESingleLayer{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    β::U
    tau::P
end

struct VIEBoundary{T,U,P} <: BoundaryOperator_ΓΩ
    gamma::T
    α::U
    tau::P
end

struct VIESingleLayer2{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    β::U
    tau::P
end

struct VIEBoundary2{T,U,P} <: BoundaryOperator_ΓΩ
    gamma::T
    α::U
    tau::P
end


struct VIEDoubleLayer{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end


struct VIEhhVolume{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end

struct VIEhhBoundary{T,U,P} <: BoundaryOperator_ΓΩ
    gamma::T
    α::U
    tau::P
end

struct VIEhhVolumegradG{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end

struct VIEhhVolumek0{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end


"""
The following operators are for test purposes only
"""

struct DVIE_TestOp1{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end

struct DVIE_TestOp2{T,U,P} <: BoundaryOperator_ΓΩ
    gamma::T
    α::U
    tau::P
end

struct DVIE_TestOp3{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end

struct EVIE_TestOp1{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end

struct EVIE_TestOp2{T,U,P} <: BoundaryOperator_ΓΩ
    gamma::T
    α::U
    tau::P
end

struct EVIE_TestOp3{T,U,P} <: VolumeOperator
    gamma::T
    α::U
    tau::P
end



#TODO: ... scalartype is also determined by the cell material...
scalartype(op::VIEOperator) = typeof(op.gamma)

export VIE




"""
Integrands for the operators of the EVIE and DVIE
"""
function (igd::Integrand{<:VIESingleLayer})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r) # Derivation after x (test variable) => "-" to get nablaG(r,r')
    
    Ty = igd.operator.tau(cartesian(y))

    αgreenTy = α * green * Ty
    βgradgreenTy = β * gradgreen * Ty


    _integrands(f,g) do fi, gi
        dot(fi.value, gi.value)*αgreenTy - dot(fi.divergence*gi.value, βgradgreenTy)
    end
end

function (igd::Integrand{<:VIESingleLayer2})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r) # Derivation after x (test variable) => "-" to get nablaG(r,r')
    
    Ty = igd.operator.tau(cartesian(y))

    αgradgreenTy = α * gradgreen * Ty

    _integrands(f,g) do fi, gi
        dot(fi.curl, cross(αgradgreenTy, gi.value))
    end
end

function (igd::Integrand{<:VIEBoundary})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r) # Derivation after x (test variable) => "-" to get nablaG(r,r')
    
    Ty = igd.operator.tau(cartesian(y))

    αgradgreenTy = α * gradgreen * Ty

    _integrands(f,g) do fi, gi
        fi.value*dot(gi.value, αgradgreenTy)
    end
end

function (igd::Integrand{<:VIEBoundary2})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r) # Derivation after x (test variable) => "-" to get nablaG(r,r')
    
    Ty = igd.operator.tau(cartesian(y))

    αgradgreenTy = α * gradgreen * Ty

    _integrands(f,g) do fi, gi
        dot(fi.value, cross(αgradgreenTy, gi.value))
    end
end

function integrand(viop::VIEDoubleLayer, kerneldata, tvals, tgeo, bvals, bgeo)
    gx = tvals[1]
    fy = bvals[1]

    gradG = kerneldata.gradgreen
    
    Ty = kerneldata.tau
    
    α = viop.α

    t = α * dot(gx, cross(gradG, Ty*fy))
end


"""
Integrands for the operators of the Lippmann-Schwinger Volume Integral Equation
"""
function (igd::Integrand{<:VIEhhVolume})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)

    Ty = igd.operator.tau(cartesian(y))

    αGTy = α * green * Ty

    _integrands(f,g) do fi, gi
        dot(fi.gradient, αGTy*gi.gradient)
    end
end

function (igd::Integrand{<:VIEhhBoundary})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)

    Ty = igd.operator.tau(cartesian(y))

    αGTy = α * green * Ty

    nx = x.patch.normals[1]

    _integrands(f,g) do fi, gi
        fi.value*dot(nx, αGTy*gi.gradient)
    end
end

function (igd::Integrand{<:VIEhhVolumek0})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)

    Ty = igd.operator.tau(cartesian(y))

    αgreenTy = α * green * Ty

    _integrands(f,g) do fi, gi
        dot(fi.value, gi.value)*αgreenTy
    end
end

function (igd::Integrand{<:VIEhhVolumegradG})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = (γ + iR) * green * (iR * r) # Derivation after y (trial variable) => "+" to get nabla'G(r,r')
    
    Ty = igd.operator.tau(cartesian(y))

    αgradgreenTy = α * gradgreen * Ty

    _integrands(f,g) do fi, gi
        fi.value*dot(gi.gradient, αgradgreenTy)
    end
end


"""
Integrands for the operators which are for test purposes only
"""
function (igd::Integrand{<:DVIE_TestOp1})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r) # Derivation after x (test variable) => "-" to get nablaG(r,r')

    Ty = igd.operator.tau(cartesian(y))

    αgradgreenTy = α * gradgreen * Ty

    _integrands(f,g) do fi, gi
        dot(fi.value, αgradgreenTy)*gi.divergence
    end
end

function (igd::Integrand{<:DVIE_TestOp2})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)

    Ty = igd.operator.tau(cartesian(y))

    αgreenTy = α * green * Ty

    _integrands(f,g) do fi, gi
        dot(fi.value, αgreenTy)*gi.divergence
    end
end

function (igd::Integrand{<:DVIE_TestOp3})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)

    Ty = igd.operator.tau(cartesian(y))

    αgreenTy = α * green * Ty

    _integrands(f,g) do fi, gi
        fi.divergence*gi.divergence*αgreenTy
    end
end

function (igd::Integrand{<:EVIE_TestOp1})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)

    Ty = igd.operator.tau(cartesian(y))

    αgreenTy = α * green * Ty

    _integrands(f,g) do fi, gi
        dot(fi.curl, gi.value)*αgreenTy
    end
end

function (igd::Integrand{<:EVIE_TestOp2})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)

    Ty = igd.operator.tau(cartesian(y))

    αgreenTy = α * green * Ty

    _integrands(f,g) do fi, gi
        dot(fi.value, gi.value)*αgreenTy    # fi.values is n̂×NedeleccBasis means ttrace(NedeleccBasis)
    end 
end

function (igd::Integrand{<:EVIE_TestOp3})(x,y,f,g)
    α = igd.operator.α
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r) # Derivation after x (test variable) => "-" to get nablaG(r,r')

    Ty = igd.operator.tau(cartesian(y))

    αgradgreenTy = α * gradgreen * Ty

    _integrands(f,g) do fi, gi
        dot(fi.value, cross(αgradgreenTy, gi.value))
    end
end







defaultquadstrat(op::VIEOperator, tfs, bfs) = SauterSchwab3DQStrat(3,3,3,3,3,3)


function quaddata(op::VIEOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::SauterSchwab3DQStrat)

    #The combinations of rules (6,7) and (5,7 are) BAAAADDDD
    # they result in many near singularity evaluations with any
    # resemblence of accuracy going down the drain! Simply don't!
    # (same for (5,7) btw...).
    t_qp = quadpoints(test_local_space,  test_charts,  (qs.outer_rule,))
    b_qp = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))

   
    sing_qp = (SauterSchwab3D._legendre(qs.sauter_schwab_1D,0,1), 
               SauterSchwab3D._shunnham2D(qs.sauter_schwab_2D),
               SauterSchwab3D._shunnham3D(qs.sauter_schwab_3D),
               SauterSchwab3D._shunnham4D(qs.sauter_schwab_4D),)


    return (tpoints=t_qp, bpoints=b_qp, sing_qp=sing_qp)
end



function _hits(τ, σ)
    T = coordtype(τ)
    hits = 0
    dtol = 1.0e3 * eps(T)
    idx_t = Int64[]
    idx_s = Int64[]
    sizehint!(idx_t,4)
    sizehint!(idx_s,4)

    for (i,t) in enumerate(τ.vertices)
        for (j,s) in enumerate(σ.vertices)
            #d2 = LinearAlgebra.norm_sqr(t-s)
            d = norm(t-s)
            #dmin2 = min(dmin2, d2)
            # if d2 < dtol
            if d < dtol
                push!(idx_t,i)
                push!(idx_s,j)
                hits +=1
                break
            end
        end
    end

    return hits, idx_t, idx_s
end




function quadrule(op::VolumeOperator, g::RefSpace, f::RefSpace, 
    i, τ::CompScienceMeshes.Simplex{<:Any, 3}, 
    j, σ::CompScienceMeshes.Simplex{<:Any, 3}, 
    qd, qs::SauterSchwab3DQStrat)

    hits, idx_t, idx_s = _hits(τ, σ)

    @assert hits <= 4

    hits == 4 && return SauterSchwab3D.CommonVolume6D_S(SauterSchwab3D.Singularity6DVolume(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[4]))
    hits == 3 && return SauterSchwab3D.CommonFace6D_S(SauterSchwab3D.Singularity6DFace(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 2 && return SauterSchwab3D.CommonEdge6D_S(SauterSchwab3D.Singularity6DEdge(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3],qd.sing_qp[4]))
    hits == 1 && return SauterSchwab3D.CommonVertex6D_S(SauterSchwab3D.Singularity6DPoint(idx_t,idx_s),qd.sing_qp[3])

    return DoubleQuadRule(qd[1][1,i], qd[2][1,j])
end



# 5D integral: ∫∫_Γ ∫∫∫_Ω
function quadrule(op::BoundaryOperator_ΓΩ, g::RefSpace, f::RefSpace, 
    i, τ::CompScienceMeshes.Simplex{<:Any, 2}, 
    j, σ::CompScienceMeshes.Simplex{<:Any, 3}, 
    qd, qs::SauterSchwab3DQStrat)

    hits, idx_t, idx_s = _hits(τ, σ)
   
    @assert hits <= 3

    hits == 3 && return SauterSchwab3D.CommonFace5D_S(SauterSchwab3D.Singularity5DFace(idx_s,idx_t),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 2 && return SauterSchwab3D.CommonEdge5D_S(SauterSchwab3D.Singularity5DEdge(idx_s,idx_t),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 1 && return SauterSchwab3D.CommonVertex5D_S(SauterSchwab3D.Singularity5DPoint(idx_s,idx_t),(qd.sing_qp[3],qd.sing_qp[2]))

    return DoubleQuadRule(qd[1][1,i], qd[2][1,j])
end


# 5D integral: ∫∫∫_Ω ∫∫_Γ
function quadrule(op::BoundaryOperator_ΩΓ, g::RefSpace, f::RefSpace, 
    i, τ::CompScienceMeshes.Simplex{<:Any, 3}, 
    j, σ::CompScienceMeshes.Simplex{<:Any, 2}, 
    qd, qs::SauterSchwab3DQStrat)

    hits, idx_t, idx_s = _hits(τ, σ)
   
    @assert hits <= 3

    hits == 3 && return SauterSchwab3D.CommonFace5D_S(SauterSchwab3D.Singularity5DFace(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 2 && return SauterSchwab3D.CommonEdge5D_S(SauterSchwab3D.Singularity5DEdge(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 1 && return SauterSchwab3D.CommonVertex5D_S(SauterSchwab3D.Singularity5DPoint(idx_t,idx_s),(qd.sing_qp[3],qd.sing_qp[2]))

    return DoubleQuadRule(qd[1][1,i], qd[2][1,j])
end
using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using SphericalScattering
using Test
using Plots


BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = BEAST.SauterSchwab3DQStrat(3,3,3,3,3,3)



r = 1.0
h = 0.5#0.25
mesh = CompScienceMeshes.tetmeshsphere(r,h)
bnd = boundary(mesh)

#Mesh
#v_list = [point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0), point(0.0, 1.0, 0.0), point(0.0, 0.0, 1.0), point(0.0, 0.0, -1.0)]
#indΩ1 = SVector(2,3,4,1) #LHR optimal
#indΩ2 = SVector(3,2,5,1)
#indΩ1 = SVector(1,2,4,3) #eigentlich erlaubt, macht aber massive Probleme bei evie bzgl Symm
#indΩ2 = SVector(1,5,2,3)
#f_list = [indΩ1, indΩ2]
#mesh=Mesh(v_list, f_list)
#bnd=boundary(mesh)


## LSVIE

X = lagrangec0d1(mesh; dirichlet = false)

@show numfunctions(X)
strc = X -> strace(X,bnd)

# VIE operators

function tau(x::SVector{U,T}) where {U,T}
    return 1.0 + x[1]^2*(x[2] - 2.0) + (x[3] - 0.5)^2
end
τ = tau

I, V, B =  Identity(), VIE.hhvolume(tau = τ, wavenumber = 0.0), VIE.hhboundary(tau = τ, wavenumber = 0.0)
Y = VIE.hhvolumegradG(tau = τ, wavenumber = 0.0)


Z_I = assemble(I, X, X)

Z_V = assemble(V, X, X)
Z_B = assemble(B, strc(X), X)

Z_Y = assemble(Y, X, X)

Z_version1 = Z_I + Z_V + Z_B
Z_version2 = Z_I + Z_Y

norm(Z_version1 - Z_version2)/norm(Z_version2)
@test norm(Z_version1 - Z_version2)/norm(Z_version2) < 0.008



#DIFF = log.(abs.(Z_version1-Z_version2) .+eps(1.0))
#yflip!(Plots.heatmap(DIFF, title = "2D Rasterplot der Differenzmatrix"); size=(1200,1000))




## Points inside Sphere !check boundary!
###########################################
###############################################
## IST DIE IMPLEMENTIERUNG VON SphericalScattering FALSCH??????????? PRUEFEN::: siehe Plot...
## Was ist mit dem TestOp3??????
###############################
###############################
#############################

ins_range_r=range(r * (1/100), stop = r * (80/100), length = 11) #!!!Anfang und Ende kritisch
ins_range_φ=range(0.0,stop=2*pi,length=31)
ins_range_θ=range(pi*(1/100),stop=pi*(99/100),length=18) # Punkte ausschließen!!!
points_inside = [point(r * cos(φ) * sin(θ), r * sin(φ) * sin(θ), r * 
cos(θ)) for r in ins_range_r for θ in ins_range_θ for φ in ins_range_φ]


ϵ0 = 1.0
ϵ1 = 4.0
μ0 = 1.0
ω = 1.0
κ0 = ω * √(ϵ0*μ0)
κ1 = ω * √(ϵ1*μ0)
#η = √(μ/ϵ)
ϵ_r = ϵ1

amplitude = 1.0
polarization = x̂
direction = ẑ

#REFERENZ SphericalScattering
sp = DielectricSphere(; radius=r, filling=Medium(ϵ1, μ0))
ex = planeWave(; embedding=Medium(ϵ0, μ0), frequency=ω/(2*pi))#Ebene Welle zur Anregung: E ist die BEAST.PlaneWaveVIE
EF₁ = field(sp, ex, ElectricField(points_inside))


## evie

function tau(x::SVector{U,T}) where {U,T}
    ϵ_r-1.0
end

ttrc = X -> ttrace(X,bnd)

X = nedelecc3d(mesh)
@show numfunctions(X)

χ = tau
I2 = Identity()
K2 = VIE.singlelayer2(wavenumber=κ0, tau=χ)
B2 = VIE.boundary2(wavenumber=κ0, tau=χ)

M_I2 = assemble(I2, X, X)
M_K2 = assemble(-K2, X, X)
M_B2 = assemble(B2, ttrc(X), X)

M_evie = ϵ_r*M_I2 + M_K2 + M_B2

E = VIE.planewave(direction=direction, polarization=polarization, wavenumber=κ0)
b = assemble(E, X)

u_evie = M_evie \ b

#Z = range(-1,1,length=100)
#Y = range(-1,1,length=100)
#nfpoints = [point(0,y,z) for  y in Y, z in Z]

Enear_evie = BEAST.grideval(points_inside,u_evie, X)

err_evie = norm(EF₁ - Enear_evie) / norm(EF₁)
@show err_evie


## dvie 

function tau(x::SVector{U,T}) where {U,T}
    1.0-1.0/ϵ_r
end

ntrc = X -> ntrace(X,bnd)

X = nedelecd3d(mesh)
@show numfunctions(X)

χ = tau
I = Identity()
K = VIE.singlelayer(wavenumber = κ0, tau = χ)
B = VIE.boundary(wavenumber = κ0, tau = χ)

M_I = assemble(I, X, X)
M_K = assemble(-K, X, X)
M_B = assemble(-B, ntrc(X), X)

M_dvie = (1.0/ϵ_r)*M_I + M_K + M_B

E = VIE.planewave(direction=ẑ, polarization=x̂, wavenumber=κ0)
b = assemble(E, X)

u_dvie = M_dvie \ b

#Z = range(-1,1,length=100)
#Y = range(-1,1,length=100)
#nfpoints = [point(0,y,z) for  y in Y, z in Z]

Enear_dvie = (1/ϵ_r).*BEAST.grideval(points_inside,u_dvie, X)

err_dvie = norm(EF₁ - Enear_dvie) / norm(EF₁)
@show err_dvie


##

norm(Enear_dvie - Enear_evie)/norm(Enear_evie)

# Gab es sicher zwei versionen??? JA 
# A) Vergleich mit SphericalScattering
# B) Term herleitung von dyade - dort könnte es alternativen geben für material = konst zum vergleichen....
# C) Dyade war ja kritisch weil es da einen versteckten teil gibt...




## dvie testop - Fraglich warum das nicht geht

function tau(x::SVector{U,T}) where {U,T}
    return 1.0 + x[1]^2*(x[2] - 2.0) + (x[3] - 0.5)^2
end

ntrc = X -> ntrace(X,bnd)

X = nedelecd3d(mesh)
@show numfunctions(X)

TestOp1 = BEAST.VIE_dvie_TestOp1(1.0, 1.0, tau)
TestOp2 = BEAST.VIE_dvie_TestOp2(1.0, 1.0, tau)
TestOp3 = BEAST.VIE_dvie_TestOp3(1.0, -1.0, tau)

M_1 = assemble(TestOp1, X, X)
M_2 = assemble(TestOp2, ntrc(X), X)
M_3 = assemble(TestOp3, X, X)
M_0 = M_2 + M_3

@show norm(M_0 - M_1)/norm(M_1) 

norm(M_1)
norm(M_2)
norm(M_3)
norm(M_0)

abs.((M_0 - M_1)./M_1)


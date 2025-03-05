using SphericalScattering
using BEAST
using CompScienceMeshes
using StaticArrays
using Plotly
using LinearAlgebra
using Plots
using PlotlyJS

# TESTE ALTE IMPLEMENTIERUNG von CEDRIC!!!!!!!!!!!!!!!!!!!

ϵ0 = 1.0
ϵ1 = 4.0
μ0 = 1.0
ω = 1.0
κ0 = ω * √(ϵ0*μ0)
#η = √(μ/ϵ)
ϵ_r = ϵ1

amplitude = 1.0
polarization = x̂
direction = ẑ



r = 1.0
h = 0.4 #0.25
mesh = CompScienceMeshes.tetmeshsphere(r,h)
bnd = boundary(mesh)


p = patch(mesh)#, fcr=nothing; caxis=nothing, showscale=true, color="red", kwargs...)
p = CompScienceMeshes.wireframe(mesh)
PlotlyJS.plot(p)


## REFERENZ SphericalScattering - Observe the potential on a x-line
sp = DielectricSphere(; radius=r, filling=Medium(ϵ1, μ0))
ex = planeWave(; amplitude=amplitude, polarization = polarization, direction = direction, embedding=Medium(ϵ0, μ0), frequency=ω/(2*pi))#Ebene Welle zur Anregung: E ist die BEAST.PlaneWaveVIE
y0 = 0.0
z0 = 0.0
x_range = range(-0.99*r, stop = 0.99*r, length = 201)
points_x = [point(x, y0, z0) for x in x_range]
x = collect(x_range)
E_on_x = field(sp, ex, ElectricField(points_x))


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

E = VIE.planewave(amplitude=amplitude, direction=direction, polarization=polarization, wavenumber=κ0)
b = assemble(E, X)

u_evie = M_evie \ b

Eevie_on_x = BEAST.grideval(points_x, u_evie, X)



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

E = VIE.planewave(amplitude=amplitude, direction=direction, polarization=polarization, wavenumber=κ0)
b = assemble(E, X)

u_dvie = M_dvie \ b

Edvie_on_x = (1/ϵ_r).*BEAST.grideval(points_x, u_dvie, X)


## Plot

Plots.plot(x, norm.(E_on_x), label = "SphericalScattering")
Plots.plot!(x, norm.(Eevie_on_x), label = "evie")
Plots.plot!(x, norm.(Edvie_on_x), label = "dvie")
xlims!(-1.0, 1.0)
#ylims!(-0.3, 0.0)
title!("|ElectricField|")
xlabel!("x")




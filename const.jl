#this file houses the constants used in the calculations

#conversion of dimensions
cm = 0.01
nm = 1e-9
μm = 1e-6

#some fundamental constants
c_0 = 299792458 # speed of light in vacuum [m/s]
h = 6.626e-34  # Planck constant [Js]

#constants for the calculation
steps_time = 1000
steps_crystal = 1000

#include here some material constants, mainly for comparison

#structs to contain data
#now make the structs using the abstract types

# how to use struct xxx -> t = xxx ()...
struct material
    #content
    Length::Float64
    Wavelength::Array{Float64}
    σ_abs::Array{Float64}
    σ_em::Array{Float64}
    τ_fluo::Float64
    N1per::Float64
    Doping::Float64
    N_gradient::Array{Float64}
end

struct laserparameters
    #parameter of the potential laser
    Wavelength::Array{Float64}
    Intensity_spectral::Array{Float64}
    Time::Array{Float64}
    Intensity_temporal::Array{Float64}
end

mutable struct ray_ABCD
    A::Complex{Float64}
    B::Complex{Float64}
    C::Complex{Float64}
    D::Complex{Float64}
end

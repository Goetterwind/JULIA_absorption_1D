#this file houses the constants used in the calculations

#conversion of dimensions
const cm = 0.01
const nm = 1e-9
const μm = 1e-6

#some fundamental constants
const c_0 = 299792458 # speed of light in vacuum [m/s]
const h = 6.626e-34  # Planck constant [Js]

#constants for the calculation
steps_time = 300
steps_crystal = 100

#spectral limits
λ_s = 900 * nm
λ_e = 1100 * nm
Δλ = 1 *nm
λ_p = 940 *nm

#wavelength Array
# λ_array = LinRange(λ_s, λ_e, 201)
λ_array = λ_s: Δλ :λ_e

#include here some material constants, mainly for comparison

#structs to contain data
#now make the structs using the abstract types

# mutable structs are slower!

# how to use struct xxx -> t = xxx ()...
struct Material
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

# generate a new struct for the individual output after manipulation
mutable struct Laserparameters
    #parameter of the potential laser
    Wavelength::Array{Float64}
    Intensity_spectral::Array{Float64}
    Time::Array{Float64}
    Intensity_temporal::Array{Float64}
end

#  speed increase by generating new rays
mutable struct Ray_ABCD
    A::Complex{Float64}
    B::Complex{Float64}
    C::Complex{Float64}
    D::Complex{Float64}
end

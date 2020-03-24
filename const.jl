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
struct material
    #content
    Length::Float32
    Wavelength::Array{Float32}
    σ_abs::Array{Float32}
    σ_em::Array{Float32}
    τ_fluo::Float32
    N1per::Float32
    Doping::Float32
    N_gradient::Array{Float32}
end

struct laserparameters
    #parameter of the potential laser
    Wavelength::Array{Float32}
    Intensity_spectral::Array{Float32}
    Time::Array{Float32}
    Intensity_temporal::Array{Float32}
end

# main file to run the project
# this is a preliminary version
#
#
# Daniel Albach
#@time begin

#using section
using ProgressMeter
using ProgressLogging
using DelimitedFiles
using Plots
using Interpolations
using BenchmarkTools


# fundamental constants and functions
include("const.jl")
include("functions.jl")

#load cross sections, otherwise use the const.jl values, which are an array as well
#loading using readdlm out of DelimitedFiles

#read the absorption data
filename = "JK_CaF300Ka.txt"
subpath = "original_code"
filepath = joinpath(@__DIR__,subpath,filename)

spectra_abs = readdlm(filepath)

# clear the NaNs that are potentially hidden inside the spectrum and generate the interpolation function
spectra_abs[:,2] = replace!(spectra_abs[:,2], NaN => 0)
spectra_abs_ip = interpolate((spectra_abs[:,1],), spectra_abs[:,2], Gridded(Linear()))

#read the fluorescence data
filename = "JK_CaF300Kf.txt"
subpath = "original_code"
filepath = joinpath(@__DIR__,subpath,filename)

spectra_flu = readdlm(filepath)

# clear the NaNs that are potentially hidden inside the spectrum and generate the interpolation function
spectra_flu[:,2] = replace!(spectra_flu[:,2], NaN => 0)
spectra_flu_ip = interpolate((spectra_flu[:,1],), spectra_flu[:,2], Gridded(Linear()))

# now make a nice plot of the spectra overlaying each other
#= pll = plot(spectra_abs[:,1],spectra_abs[:,2], ylims=(0,1e-20))
plot!(spectra_flu[:,1],spectra_flu[:,2], ylims=(0,1e-20)) =#

# display(pll)

# here the actual code starts
# for sure some monochromatic one at first

#println("constant speed of light $c_light")
#@showprogress 1 "Computing..."
#scope is a bit of a pain - use let end or define as globals
#let
#=     #ProgreesMeter call
    p=Progress(steps_crystal, dt=0.05, color=:grey)

    a=10
    b=20

    #= @progress for i in 1:steps_crystal
        # sleep(0.005)
        a,b = test_function(i,steps_crystal)
        ProgressMeter.next!(p)
        # sleep(0.005)
    end =#

    println(a);
    println(b);
    
    p,a,b = iter(p,a,b)
    println(a);
    println(b); =#
#end

# so do the monochromatic version first - generate a vector with zeroes for the beta-distribution
β_vec = zeros(1,steps_crystal);
β_inter = zeros(1,steps_crystal-1);

# define the cross sections for now
σ_ap = 0.7e-20; #cm^-2
σ_ep = 0.22e-20; #cm^-2

# material size etc.
crys_l = 0.7; #cm
crys_d = 2 * 1.388e20;

crys_step = crys_l /(steps_crystal-1);

# now do a constant pump intensity for the sake of simplicity, which will be the later time dependent
I_pump = 16e3; # W/cm^2
τ_fluo =0.94e-3; # s

pump_dur = 1e-3;
δt = pump_dur / (steps_time-1);

# this is the initial pump intensity vector along the crystal axis
pump_vec = ones(1,steps_crystal);

for itime in 1:steps_time
    pump_vec[1] = I_pump;

    # later add the pump recycling and the multipump version, for now a simple onesided, no gradient of the doping yet
    for icrys in 1:steps_crystal-1
        # in order to estimate the absorption factor, we have the classic fence problem
        # interpolate fromt the former the initial distribution the new β ;)
        local β_inter = (β_vec[icrys] + β_vec[icrys+1])/2;
        pump_vec[icrys+1] = pump_vec[icrys] * exp( -(σ_ap-β_inter*(σ_ap+σ_ep))*crys_d*crys_step);
    end

    # display(plot(pump_vec[1,:]))

    # now we can integrate to get the new β distribution using the differential equation, explicit solution of the diffeq
    local A1vec = σ_ap.*pump_vec./(h*c_0/λ_p);
    local C1vec = (σ_ap+σ_ep).*pump_vec./(h*c_0/λ_p).+1/τ_fluo;
    global β_vec = A1vec./C1vec .* (1 .-exp.(-C1vec.*δt)) .+ β_vec.*exp.(-C1vec.*δt);

    # the very same can be done in a relatively longer version
    #= for ibeta in 1:steps_crystal
        # A1 = σ_ap*pump_vec[ibeta]/(h*c_0/λ_p);
        # C1 = (σ_ap+σ_ep)*pump_vec[ibeta]/(h*c_0/λ_p)+1/τ_fluo;
        A1 = A1vec[ibeta];
        C1 = C1vec[ibeta]
        β_vec[ibeta] = A1/C1 * (1-exp(-C1*δt)) + β_vec[ibeta]*exp(-C1*δt);
    end =#

    # display(plot(β_vec[1,:]))
end

display(plot(β_vec[1,:]))
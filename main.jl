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


# fundamental constants
include("const.jl")
# necessary functions
include("functions.jl")

#load cross sections, otherwise use the const.jl values, which are an array as well
#loading using readdlm out of DelimitedFiles

#read the absorption data
filename = "JK_CaF300Ka.txt"
subpath = "original_code"
filepath = joinpath(@__DIR__,subpath,filename)

spectra_abs = readdlm(filepath)

# clear the NaNs that are potentially hidden inside the spectrum
spectra_abs[:,2] = replace!(spectra_abs[:,2], NaN => 0)
spectra_abs_ip = interpolate((spectra_abs[:,1],), spectra_abs[:,2], Gridded(Linear()))

#read the fluorescence data
filename = "JK_CaF300Kf.txt"
subpath = "original_code"
filepath = joinpath(@__DIR__,subpath,filename)

spectra_flu = readdlm(filepath)

# clear the NaNs that are potentially hidden inside the spectrum
spectra_flu[:,2] = replace!(spectra_flu[:,2], NaN => 0)

# do some linear interpolation in the interesting range
spectra_flu_ip = interpolate((spectra_flu[:,1],), spectra_flu[:,2], Gridded(Linear()))

# now make a nice plot of the spectra overlaying each other
plot(spectra_abs[:,1],spectra_abs[:,2], ylims=(0,1e-20))s
plot!(spectra_flu[:,1],spectra_flu[:,2], ylims=(0,1e-20))

#println("constant speed of light $c_light")
#@showprogress 1 "Computing..."
#scope is a bit of a pain - use let end or define as globals
#let
    #ProgreesMeter call
    # p=Progress(steps_crystal, dt=0.05, color=:grey)

    a=0
    b=0

    @progress for i in 1:steps_crystal
        #sleep(0.1)
        global a
        global b
        a,b = test_function(i,steps_crystal)
        # ProgressMeter.next!(p)
        sleep(0.02)
    end
    println(a);
#end

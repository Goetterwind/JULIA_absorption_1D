# main file to run the project
# this is a preliminary version
#
#
# Daniel Albach
#@time begin

#using section
using ProgressMeter
using DelimitedFiles
using Plots

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

# spectra_abs = readdlm(filepath)

#read the fluorescence data
filename = "JK_CaF300Kf.txt"
subpath = "original_code"
filepath = joinpath(@__DIR__,subpath,filename)

# spectra_flu = readdlm(filepath)

#println("constant speed of light $c_light")
#@showprogress 1 "Computing..."
#scope is a bit of a pain - use let end or define as globals
#let
    #ProgreesMeter call
    p=Progress(steps_crystal, dt=0.05, color=:grey)

    a=0
    b=0

    for i in 1:steps_crystal
        #sleep(0.1)
        global a
        global b
        a,b = test_function(i,steps_crystal)
        ProgressMeter.next!(p)
        sleep(0.1)
    end
    println(a);
#end

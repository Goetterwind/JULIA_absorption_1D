# main file to run the project
# this is a preliminary version
#
#
# Daniel Albach
#@time begin

#using section
using ProgressMeter

# fundamental constants
include("const.jl")
# necessary functions
include("functions.jl")

#println("constant speed of light $c_light")
#@showprogress 1 "Computing..."
#scope is a bit of a pain - use let end or define as globals
#let
    #ProgreesMeter call
    p=Progress(steps_crystal, dt=1.0, color=:grey)

    a=0
    b=0

    for i in 1:steps_crystal
        #sleep(0.1)
        global a
        global b
        a,b = test_function(i,steps_crystal)
        ProgressMeter.next!(p)
    end
    println(a);
#end

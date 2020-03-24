# main file to run the project
# this is a preliminary version
#
#
# Daniel Albach

#using section
using ProgressMeter

# fundamental constants
include("const.jl")

# necessary functions
include("functions.jl")

#println("constant speed of light $c_light")
@showprogress 1 "Computing..." for i in 1:50
    sleep(0.1)
end
println("Hello")

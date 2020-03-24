# main file to run the project
# this is a preliminary version
#
#
# Daniel Albach
#@time begin

#using section
using ProgressMeter
p=Progress(steps_crystal, dt=1.0, color=:grey)
# fundamental constants
include("const.jl")

# necessary functions
include("functions.jl")

#println("constant speed of light $c_light")
#@showprogress 1 "Computing..."
for i in 1:steps_crystal
    sleep(0.01)
    ProgressMeter.next!(p)
end

#end

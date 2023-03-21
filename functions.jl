#this is the accumulation of the different functions needed for the main file

#test function in order to elarn how to handle functions in Julia
#this function returns secveral "variables", but in reality they are tuples
function test_function(a,b)
    return a+b, a*b
end

# this function does a reverse ray trace to get the initial distribution of the pump field
# this has to be meant as a free space ray trace propagation
# input is a (sum of) phase distribution and an intensity distribution
# input is number of desired rays
# output are the rays as vectors
# in case of an incoherent source MLA, it is just the MLA property and Fourier lens
function reverse_3D_tracing(args)
    # body
end

# Gaussian ABCD propagation including the new distribution and the phase of the center beam
# a good test would be to get a Poissons spot example for diffraction
function gaussian_prop_ABCD(args)
    # body
end

# this function generates a ray (either geometric or Gaussian)
function generate_ray(args)
    # body
end

function β_int(args)
    #this is a copy of the old Matlab version of the estimation of β

end

function iter(p,a,b)
    @progress for i in 1:steps_crystal
        # sleep(0.005)
        a,b = test_function(i,steps_crystal)
        ProgressMeter.next!(p)
    end
    return p,a,b
end

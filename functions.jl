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

function get_pump_vec(pump_vec,β_vec,p_dir)
    if p_dir == 1
        #do the pos
        pump_vec = do_pump_vec(pump_vec,β_vec);
        pump_ret = pump_vec[steps_crystal];
    elseif p_dir == -1
        #do the neg
        # pump_vec = reverse(pump_vec);
        β_vec = reverse(β_vec);
        pump_vec = do_pump_vec(pump_vec,β_vec);
        pump_ret = pump_vec[steps_crystal];
        pump_vec = reverse(pump_vec);
    else
        #this is 0 case
        # do nothing and return the '0-vector'
        pump_vec = zeros(1,steps_crystal);
        pump_ret = 0;
    end
    # println(pump_ret);
    return pump_vec, pump_ret
end

function do_pump_vec(pump_vec,β_vec)
    for icrys in 1:steps_crystal-1
        # in order to estimate the absorption factor, we have the classic fence problem
        # interpolate fromt the former the initial distribution the new β ;)
        # crys_d later has to be changed to a vector for the doping concentration
        local β_inter = (β_vec[icrys] + β_vec[icrys+1])/2;
        pump_vec[icrys+1] = pump_vec[icrys] * exp( -(σ_ap-β_inter*(σ_ap+σ_ep))*crys_d*crys_step);
    end
    return pump_vec
end

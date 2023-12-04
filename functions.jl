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

function beta_int(s_ap::Float64, s_ep::Float64, l_p::Float64, beta_v::Matrix{Float64}, pump_v::Matrix{Float64}, dt::Float64)
    #this is a copy of the old Matlab version of the estimation of β
    β=beta_v;
    A1vec = s_ap.*pump_v./(h*c_0/l_p);
    C1vec = (s_ap+s_ep).*pump_v./(h*c_0/l_p).+1/τ_fluo;
    β = A1vec./C1vec .* (1 .-exp.(-C1vec.*dt)) .+ β.*exp.(-C1vec.*dt);
    return β
end

function iter(p,a,b)
    @progress for i in 1:steps_crystal
        # sleep(0.005)
        a,b = test_function(i,steps_crystal)
        ProgressMeter.next!(p)
    end
    return p,a,b
end

function get_pump_vec(pv,bv,pdir::Int64)
    if pdir == 1
        #do the pos
        pv = do_pump_vec(pv,bv);
        pump_ret = pv[steps_crystal];
    elseif pdir == -1
        #do the neg
        # pump_vec = reverse(pump_vec);
        bv = reverse(bv);
        pv = do_pump_vec(pv,bv);
        pump_ret = pv[steps_crystal];
        pv = reverse(pv);
    else
        #this is 0 case
        # do nothing and return the '0-vector'
        pv = zeros(steps_crystal);
        pump_ret = 0;
    end
    # println(pump_ret);
    return pv, pump_ret
end

function do_pump_vec(pump_vec,β_vec)
    for icrys in 1:steps_crystal-1
        # in order to estimate the absorption factor, we have the classic fence problem
        # interpolate fromt the former the initial distribution the new β ;)
        # crys_d later has to be changed to a vector for the doping concentration
        β_inter = (β_vec[icrys] + β_vec[icrys+1])/2;
        pump_vec[icrys+1] = pump_vec[icrys] * exp( -(σ_ap-β_inter*(σ_ap+σ_ep))*crys_d*crys_step);
    end
    return pump_vec
end

#= function pump_next(pump_vec, cind, I_pv, pump_ret)
    if I_pv == 1
        pump_vec[cind+1,1]=pump_ret;
    elseif I_pv == -1
        pump_vec[cind+1,1]=pump_ret;
    elseif I_pv == 0
        pump_vec[cind+1,1]=0;
    end
    return pump_vec
end =#

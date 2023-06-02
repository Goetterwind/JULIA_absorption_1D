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
#using BenchmarkTools


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

# so do the monochromatic version first - generate a vector with zeroes for the beta-distribution
# it is not good practice to use a matrix as an array -> initialize as an array
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
I_pump = [10e3 10e3]; # W/cm^2
τ_fluo =0.94e-3; # s

pump_dur = 1e-3;
δt = pump_dur / (steps_time-1);

# how to define the brm and the multiple pumps?
# directions are indicated by the sign [], the '0' would ignore it in case

I_pv = [-1 1; 1 -1] #brm basically

# this is the initial pump intensity vector along the crystal axis
pump_vec = zeros(size(I_pv,1)*size(I_pv,2),steps_crystal);


# some tip for debugging: if steps_time ~= 1 then you will see saturation effects!

# how to make this braod band? We basically have to add the spectral distribution beforehand and give it to the next
# iteration instead of a single intensity value
# better seperate the monochromatic from the polychromatic version!

# What is still missing? We still need to itroduce the doping vector and reflectivities after each pass!

for itime in 1:steps_time
    # later add the pump recycling and the multipump version, for now a simple onesided, no gradient of the doping yet
    # we now have to go through the size of the I_pv in its two dimensions
    # flip the arrays using 'reverse'?

    # size(I_pv,1) gives the amount of available pumps
    # size(I_pv,2) gives the amount of rountrips
    # you have to initialize the pump vector with 0's and give the first 

    # but this can be turned into a function call
    for ip in 1:size(I_pv,1)
        cind = size(I_pv,2)*(ip-1)+1;
        pump_vec[cind,1] = I_pump[ip];
        for ir in 1:size(I_pv,2)
            cind = size(I_pv,2)*(ip-1)+ir;
            # make a local variable or put this whole lopp straight into a function
            global pump_vec[cind,:], pump_ret = get_pump_vec(pump_vec[cind,:],β_vec,I_pv[ip,ir])
            if ir == size(I_pv,2)
                continue
            end
            # turning this into a function does not give any advantage
             if I_pv[ip,ir] == 1
                pump_vec[cind+1,1]=pump_ret;
            elseif I_pv[ip,ir] == -1
                pump_vec[cind+1,1]=pump_ret;
            elseif I_pv[ip,ir] == 0
                pump_vec[cind+1,1]=0;
            end
        end
    end

    pump_v = sum(pump_vec, dims=1);

    # display(plot(pump_vec[1,:]))

    # now we can integrate to get the new β distribution using the differential equation, explicit solution of the diffeq
    A1vec = σ_ap.*pump_v./(h*c_0/λ_p);
    C1vec = (σ_ap+σ_ep).*pump_v./(h*c_0/λ_p).+1/τ_fluo;
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

display(plot(β_vec[1,:],ylimits=(0, 0.5)))

#= fig = plot()
for i in 1:size(pump_vec,1)
    plot!(pump_vec[i,:], ylimits=(0, 16e3));
end
display(fig) =#
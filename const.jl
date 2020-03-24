#this file houses the constants used in the calculations

#conversion of dimensions
cm = 0.01
nm = 1e-9
μm = 1e-6

#some fundamental constants
c_0 = 299792458 # speed of light in vacuum [m/s]
h = 6.626e-34  # Planck constant [Js]

#include here some material constants, mainly for comparison

#Yb:YAG
N1per = 1.388e20

#structs to contain data

struct crystal
    #content
    length
end

struct laserparameters
    #content of the different lasers for pump and extraction
    wavelength
end

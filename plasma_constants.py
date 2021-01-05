
######## Constant defined here ########
pi = 3.1415926535897932384626
q0 = 1.602176565e-19 # C
m0 = 9.10938291e-31 # kg
v0 = 2.99792458e8 # m/s^2
kb = 1.3806488e-23 # J/K
mu0 = 4.0e-7*pi # N/A^2
epsilon0 = 8.8541878176203899e-12  # F/m
wavelength = 0.8e-6
frequency = v0*2*pi/wavelength
# Normalized units
exunit = m0*v0*frequency/q0
bxunit = m0*frequency/q0 
denunit = frequency**2*epsilon0*m0/q0**2
curunit = q0*denunit*v0
eneunit = m0*v0*v0
momunit = m0*v0
temunit = kb/q0 
# prefix
giga = 1e9
mega = 1e6
micron = 1e-6
nano = 1e-9
femto = 1e-15


# =============================================================================
# Post-MS Solar System Orbital Evolution (REBOUND)

# Authors: Noah Ferich, Stanley A. Baronett
# Last Updated: 7/16/2019
# Version: ***
# =============================================================================
import mesa_reader as mr
import time
import rebound
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# =============================================================================
# Import MESA history.data using PyMesaReader
# =============================================================================
# Make MesaData objects from a history files
h = mr.MesaData('Documents/Programs/1M_pre_ms_to_wd_copy/LOGS/history_to_end_agb.data')
h1 = mr.MesaData('Documents/Programs/1M_pre_ms_to_wd_copy/LOGS/history_to_wd.data')

# Input sequence of data files into a, m, and logR arrays
# for the age (in years), mass (in Msun units) and radius (in log10 Rsun units)

# SOLAR AGE DATA (measured in years)
# (First real age entry at approx. 1.2275e+10 yrs;
# initial 3 points used for smooth transition of 10 millions years)
a = np.array([1.226500e+10, 1.226525e+10, 1.226550e+10])
a1 = h.star_age
a2 = h1.star_age

# SOLAR MASS DATA (measured in Msun units)
# (First 3 data points used for smooth transition from 1 Msun)
m = np.array([1.0, 1.0, 1.0]) 
m1 = h.star_mass
m2 = h1.star_mass

# SOLAR RADIUS DATA (measured in log10 Rsun units)
# (First 3 data points used for smooth transition from 1 Rsun)
r = np.array([0.0, 0.0, 0.0]) 
r1 = h.log_R
r2 = h1.log_R

# Concatenate data into single arrays
a = np.append(a, [a1])
m = np.append(m, [m1])
r = np.append(r, [r1])
age = np.append(a, [a2])
mass = np.append(m, [m2])
logR = np.append(r, [r2])

# Convert logR values into AU (REBOUND default units) and store in new array called radius.
radius = np.zeros(logR.size)
for i,lrad in enumerate(logR):
    radius[i] = (10**lrad) * 0.00465047 # convert from log10(Rsun) to AU

#Creates spline functions with default boundary condition type
csm = CubicSpline(age, mass) # mass function  
csr = CubicSpline(age, radius) #radius function

#Plots data along with the mass spline
fig1, ax1 = plt.subplots()
ax1.scatter(age, mass)
ax1.scatter(np.linspace(age[0], age[len(age)-1], num = 1000), csm(np.linspace(age[0], age[len(age)-1], num = 1000)), label='Default')
ax1.legend(loc='best', ncol=1)
ax1.set_xlabel('Age of Star(Ten Billion Years)')
ax1.set_ylabel("Mass of Star (Solar Mass))")
ax1.set_title("Change in Mass of Sun-like star- Post Main Sequence to WD")



#Simulation begins here
sim = rebound.Simulation()

#Time step for the simulation- 1/20 the orbital period of mercury
sim.dt = .012045

sim.integrator = "mercurius" #changes which integrator the simulation uses

#Changes simulation and G to units of solar masses, years, and AU
sim.units = ('yr', 'AU', 'Msun') 

#Adds Sun and outer planets to simulation
sim.add("Sun") 
sim.add("Mercury")
sim.add("Venus")
sim.add("Earth")
sim.add("Mars")
sim.add("Jupiter")
sim.add("Saturn")
sim.add("Uranus")
sim.add("Neptune")

#Moves all particles to center of momentum frame
sim.move_to_com()

#Gives orbital information before the simulation begins
for orbit in sim.calculate_orbits():
    print(orbit)


time_steps = [] #list used for testing the processing time needed for each simulation step


# Function for finding the distance to a planet's pericenter
def Pericenter(a, e):
    return a * (1-e)
    


start_all = time.perf_counter() #Start of timer

#Main simulation loop
while sim.t < 30000: #determines how many years the simulation goes for
        start = time.perf_counter() #beginning of time step timer
        sim.step() #moves simulation forward a time step
        sim.particles[0].m = csm(sim.t + age[0]) #changes the mass of the Sun
        if (sim.t/.012045) % (1000/.012045) <= 1: # goes to radius loop every thousand years in the simulation
            if Pericenter(sim.particles[1].a, sim.particles[1].e) <= csr(sim.t + age[0]): # checks to see if the planet's pericenter is within the radius of the star
                sim.particles[0].m += sim.particles[1].m 
                sim.remove(1)
                print("Closest Planet Engulfed")
                print(sim.t)
        sim.integrator_synchronize() #synchronizes all changes in the particles during the simulation
        end = time.perf_counter() # end of time step timer
        time_steps.append(end - start) #appends time step time to a list 

#End of simulation timer and prints time in seconds
end_all = time.perf_counter()
print(end_all - start_all) 

print(sim.particles[0].m) #final mass of the star       


#Gives information on the orbits after the simulation
for orbit in sim.calculate_orbits():
    print(orbit)

#Plots the orbits of the planets
fig2 = rebound.OrbitPlot(sim, trails=True, unitlabel ="[AU]")


#Changes time steps from seconds to microseconds and appends them to a new list
time_steps_microsec = []
for item in time_steps:
    time_steps_microsec.append(item * 1000000)
    
#Plots the time step processing times
fig3, ax3 = plt.subplots()
plt.axis([-5, len(time_steps), 0, 1000])
ax3.set_xlabel('Step Number')
ax3.set_ylabel("Processing Time (microseconds)")
ax3.set_title("Step Processing Time- Mercurius- 10000 steps- Radius Function 1000 Years")
ax3.scatter(np.linspace(0, len(time_steps), len(time_steps)), time_steps_microsec)
plt.show()

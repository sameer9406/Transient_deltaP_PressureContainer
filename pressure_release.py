import numpy as np
import matplotlib.pyplot as plt
# Model Parameters

#Thermal Properties
P_amb = 0.995 * 10**5 #Pa
R_air = 287.5#J/kg-k
T_amb = 55+273.15#k 
rho_amb = P_amb/R_air/T_amb
P_boundary = 15 *10**5 #Pa
T_bound = (P_boundary /P_amb)**(1- 1/gamma) * T_amb
rho_boundary =  P_boundary / (R_air * T_amb)

#Input Parameters
#Inlet Characteristics
Cp_inlet= -2
A_inlet= 50 * 10**-6 * 8#m^2
Ps_loc = Cp_inlet* gamma/2 * P_amb * Mach**2 + P_amb
mdot_in = 0.7  #kg/s
V=20 #m^3 volume of the POD
A_out = 20 * 10**-6#m^2
T_bool = 1# 1 or 0


# Simulation Parameters
t = 50 #s End of Simulation Time
delta_t = 0.01 # Sampling Time
N=int(t/delta_t)# Simulation length

     #Initilisation
rho_Pod= np.zeros(N+2) # Initialization the rho  vector
dP = np.zeros(N+2)
mdot_out = np.zeros(N+2)

k=0
rho_Pod[k] = rho_amb # Initial Vaue kg/m^3
dP[k] = 0.00000001 #Pa

mdot_out[k]= 0.00000001 #Kg/s

T= np.zeros(N+2)
T[k]= T_amb#((mdot_in * T_amb + 0.15 * 373) / (mdot_in + 0.15)) 
 
P_Pod = np.zeros(N+2)
P_Pod[k]=P_amb

# Simulation
while rho_Pod[k] < rho_boundary:
    #rho_Pod[k+1] = rho_Pod[k] + delta_t * (mdot_in - mdot_out)/V
    rho_Pod[k+1] = rho_Pod[k] + delta_t * (mdot_in - mdot_out[k])/V # Continuity Equation
    P_Pod[k+1] = rho_Pod[k+1]  * R_air*T[k]                         # Ideal gas Equation
    dP[k+1] =  P_Pod[k+1] - Ps_loc + dP_outlet                      # Difference in Pressure dP
    v_out= np.sqrt(2*dP[k+1]/rho_Pod[k+1])                          # Bernoulli equation
    mdot_out[k+1] =rho_Pod[k+1] * v_out *A_out                    # mass flow equation
    
    if T_bool == 1:
        T[k+1] = (P_Pod[k+1] /P_Pod[k])**(1- 1/gamma) * T[k]        # Adiabatic Temperature Equation
    else:
        T[k+1] =T_amb
    if k == N:
       print("Iteration Complete")
       break;
    last_value = rho_Pod[k]
    last_value_dP = dP[k] 
    last_value_mdot_out = mdot_out[k]
    last_value_T = T[k]
    last_value_P_Pod = P_Pod[k]
   
    k+=1
    
# Conditioning Array
for k in range(N+2):
    if rho_Pod[k] == 0:
       rho_Pod[k] = last_value
    if dP[k] == 0:
       dP[k] = last_value_dP
    if mdot_out[k] == 0:
       mdot_out[k] = last_value_mdot_out
   
        
    if P_Pod[k] == 0:
       P_Pod[k] = last_value_P_Pod
    if T[k] == 0:
       T[k] = last_value_T

#==============================================================================
# # Plot the Simulation Results
 #Create the Time Series

#rho_Pod = rho_Pod/V
#dP = (rho_Pod*287.05*293*10**-5) - 1
#mdot_out =  mdot_out * 0.002*1000
#dP = dP*287.05*293*10**-5
#fig = plt.figure()

t_range=np.arange(0,t+2*delta_t,delta_t)

#plt.plot(t_range, dP*1000 )
plt.subplot(1, 2, 1)
if T_bool == 1:
    plt.plot(t_range,  dP/100, label='adb. Temp') 
else:
    plt.plot(t_range,  dP/100, label='const. Temp')
    plt.grid()
#plt.plot(t_range,  rho_Pod* R_air*T_amb/100)
plt.title('              || delta P ||         Area out: %d sq.cm'  %( A_out * 10**4 ))
plt.xlabel('t [s]')
plt.ylabel('dP [mbar]')
plt.legend()
plt.grid()

if T_bool == 1:
    plt.subplot(1, 2, 2)
    plt.plot(t_range,  T-T_amb)
    plt.title('|| delta T || ' )
    plt.xlabel('t [s]')
    plt.ylabel('t [°C]')


plt.grid()



#xmin = 0
#xmax = t
#ymin = 0
#ymax = 30
#plt.axis([xmin, xmax, ymin, ymax])
plt.show()
#==============================================================================

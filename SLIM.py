#!/usr/bin/env python


# Python version of SLIM.m

## Variable names changed (for consistency with thesis paper)
# 'Vcrated' changed to 'Vr' - Rated rotor velocity (m/s)
# 'pw' changed to 'Np' - No. of parallel wires
# 'wire_d' changed to 'Dw' - Diameter of selected copper wire (mm)

import numpy as np
import matplotlib as plt


## Assign Design parameters

designno = 4   # Set Design Case no.


    
if designno == 0: # Keith's Design Parameters - case no. 0

	# ElectroMagnetic constants
#        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
#        mu = 4000*mu0       # Magnetic Permeability of laminated steel    **GOOD**
        mu = 8*10**-4        # Magnetic Permeability for carbon steel    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.0105          # Rotor outer thickness (m)
        m = 3               # Number of phases	**GOOD**
        Vline = 75          # RMS line-to-line voltage (V)
        f = 250             # Supply frequency (Hz)
#         f = 62             # Supply frequency (Hz)
        p = 6               # Number of poles
        q1 = 2              # Number of slots per pole per phase	**GOOD**
        Ws = 0.038          # Width of the stator (m)	**GOOD**
        gm = 0.004          # Physical air gap (m)
        
        Srated = 0.10       # Rated slip
        Fsprime = 800       # Target thrust (N)
        Vr = 120            # Rated rotor velocity (m/s)
#         Vr = 30            # Rated rotor velocity (m/s)
        
elif designno == 1: # Keith's Design Parameters - case no. 1
        
        # ElectroMagnetic constants
        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
        mu = 4000*mu0       # Magnetic Permeability of laminated steel    **GOOD**
#         mu = 8*10**-4        # Magnetic Permeability for carbon steel    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.0105          # Rotor outer thickness (m)
        m = 3               # Number of phases	**GOOD**
        Vline = 75          # RMS line-to-line voltage (V)
        f = 250              # Supply frequency (Hz)
        p = 4               # Number of poles
        q1 = 2              # Number of slots per pole per phase	**GOOD**
        Ws = 0.038          # Width of the stator (m)	**GOOD**
        gm = 0.004           # Physical air gap (m)
        
        Srated = 0.10       # Rated slip
        Fsprime = 800       # Target thrust (N)
        Vr = 120            # Rated rotor velocity (m/s)
        
		
elif designno == 2: # Keith's Design Parameters - case no. 2
        
        # ElectroMagnetic constants
        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
        mu = 4000*mu0       # Magnetic Permeability of laminated steel    **GOOD**
#         mu = 8*10**-4        # Magnetic Permeability for carbon steel    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.0105          # Rotor outer thickness (m)
        m = 3               # Number of phases	**GOOD**
        Vline = 75          # RMS line-to-line voltage (V)
        f = 250             # Supply frequency (Hz)
        p = 2               # Number of poles
        q1 = 2              # Number of slots per pole per phase	**GOOD**
        Ws = 0.038          # Width of the stator (m)	**GOOD**
        gm = 0.004          # Physical air gap (m)
        
        Srated = 0.10       # Rated slip
        Fsprime = 800       # Target thrust (N)
        Vr = 120            # Rated rotor velocity (m/s)
        
elif designno == 3: # Keith's Design Parameters - case no. 3
        
        # ElectroMagnetic constants
        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
        mu = 4000*mu0       # Magnetic Permeability of laminated steel    **GOOD**
#         mu = 8*10**-4        # Magnetic Permeability for carbon steel    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.021           # Rotor outer thickness (m)
        m = 3               # Number of phases	**GOOD**
        Vline = 75          # RMS line-to-line voltage (V)
        f = 250             # Supply frequency (Hz)
        p = 2               # Number of poles
        q1 = 2              # Number of slots per pole per phase	**GOOD**
        Ws = 0.038          # Width of the stator (m)	**GOOD**
        gm = 0.004          # Physical air gap (m)
        
        Srated = 0.10       # Rated slip
        Fsprime = 800       # Target thrust (N)
        Vr = 120            # Rated rotor velocity (m/s)
        
elif designno == 4: # Keith's Design Parameters - case no. 4 for a puller on the track with wheels for air-gap spacing
        
        # ElectroMagnetic constants
        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.0105          # Rotor outer thickness (m)
        m = 3               # Number of phases    **GOOD**
        Vline = 150         # RMS line-to-line voltage (V)
        f = 60              # Supply frequency (Hz)
        p = 4               # Number of poles
        q1 = 2              # Number of slots per pole per phase    **GOOD**
        Ws = 0.051          # Width of the stator (m)    **GOOD**
        gm = 0.001          # Physical air gap (m)
        
        Srated = 0.10       # Rated slip
        Fsprime = 1600       # Target thrust (N)
        Vr = 120            # Rated rotor velocity (m/s)
        
elif designno == 5: # Keith's Design Parameters - case no. 4 for a puller on the track with wheels for air-gap spacing
        
        # ElectroMagnetic constants
        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.0105          # Rotor outer thickness (m)
        m = 3               # Number of phases    **GOOD**
        Vline = 150         # RMS line-to-line voltage (V)
        f = 60              # Supply frequency (Hz)
        p = 4               # Number of poles
        q1 = 2              # Number of slots per pole per phase    **GOOD**
        Ws = 0.051          # Width of the stator (m)    **GOOD**
        gm = 0.001          # Physical air gap (m)
        
        Srated = 0.05       # Rated slip
        Fsprime = 1600       # Target thrust (N)
        Vr = 120            # Rated rotor velocity (m/s)
        
elif designno == 42: # Table 4-2 Design Parameters
        
        # ElectroMagnetic constants
        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
        mu = 4000*mu0       # Magnetic Permeability of laminated steel    **GOOD**
#         mu = 8*10**-4        # Magnetic Permeability for carbon steel    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.003           # Aluminum Rotor outer thickness (m)
        m = 3               # Number of phases	**GOOD**
        Vline = 480         # RMS line-to-line voltage (V)
        f = 60              # Supply frequency (Hz)
        p = 4               # Number of poles
        q1 = 1              # Number of slots per pole per phase	**GOOD**
        Ws = 3.14           # Width of the stator (m)	**GOOD**
        gm = 0.01           # Physical air gap (m)
        
        Srated = 0.10       # Rated slip
#         Srated = 0.05       # Rated slip
        Fsprime = 8161      # Target thrust (N)
#         Fsprime = 8177      # Target thrust (N)
        Vr = 15.5           # Rated rotor velocity (m/s)
        
elif designno == 44: # Table 4-4 Design Parameters
        
        # ElectroMagnetic constants
        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
        mu = 4000*mu0       # Magnetic Permeability of laminated steel    **GOOD**
#         mu = 8*10**-4        # Magnetic Permeability for carbon steel    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.003           # Aluminum Rotor outer thickness (m)
        m = 3               # Number of phases	**GOOD**
        Vline = 480         # RMS line-to-line voltage (V)
        f = 60              # Supply frequency (Hz)
        p = 4               # Number of poles
        q1 = 1              # Number of slots per pole per phase	**GOOD**
        Ws = 3.1416         # Width of the stator (m)	**GOOD**
        gm = 0.01           # Physical air gap (m)
        
        Srated = 0.10       # Rated slip
        Fsprime = 8161      # Target thrust (N)
        Vr = 15.5           # Rated rotor velocity (m/s) 
        
elif designno == 45: # Table 4-5 Design Parameters
        
        # ElectroMagnetic constants
        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
        mu = 4000*mu0       # Magnetic Permeability of laminated steel    **GOOD**
#         mu = 8*10**-4        # Magnetic Permeability for carbon steel    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.003           # Aluminum Rotor outer thickness (m)
        m = 3               # Number of phases	**GOOD**
        Vline = 480         # RMS line-to-line voltage (V)
        f = 60              # Supply frequency (Hz)
        p = 4               # Number of poles
        q1 = 1              # Number of slots per pole per phase	**GOOD**
        Ws = 3.1416         # Width of the stator (m)	**GOOD**
        gm = 0.01           # Physical air gap (m)
        
        Srated = 0.10       # Rated slip
        Fsprime = 8171      # Target thrust (N)
        Vr = 15.5           # Rated rotor velocity (m/s) 
        
elif designno == 54: # Table 5-4 Design Parameters
        
        # ElectroMagnetic constants
        mu0 = 4*np.pi*10**-7    # Permeability of free-space    **GOOD**
        mu = 4000*mu0       # Magnetic Permeability of laminated steel    **GOOD**
#         mu = 8*10**-4        # Magnetic Permeability for carbon steel    **GOOD**
        rhow = 19.27*10**-9  # Copper volume resistivity
        rhor = 28.85*10**-9  # Capsule conductor volume resistivity
        btmax = 1.6         # Maximum allowable flux density in tooth (T)
        bymax = 1.3         # Maximum allowable flux density in yoke (T)
        J1 = 6*10**6         # Stator current density (A/m**2)
        
        # Design parameters
        d = 0.003           # Aluminum Rotor outer thickness (m)
        m = 3               # Number of phases	**GOOD**
        Vline = 480         # RMS line-to-line voltage (V)
        f = 60              # Supply frequency (Hz)
        p = 4               # Number of poles
        q1 = 1              # Number of slots per pole per phase	**GOOD**
        Ws = 3.1416         # Width of the stator (m)	**GOOD**
        gm = 0.01           # Physical air gap (m)
        
        Srated = 0.10       # Rated slip
        Fsprime = 8161      # Target thrust (N)
#         Srated = 0.05       # Rated slip
#         Fsprime = 8177      # Target thrust (N)
        Vr = 15.5           # Rated rotor velocity (m/s) 
        
#end

## Optimize for design case and calculate output performance

# Data from the PCP design procedure
V1 = Vline/np.sqrt(3)         # Rated primary RMS - Eqn 4.16
Vs = Vr/(1 - Srated)       # Sychronous velocity (m/s)	**GOOD**
tau = Vs/(2*f)             # Pole np.pitch	**GOOD**
Lambda = tau/(m*q1)        # Slot np.pitch    **GOOD**
Ls = p*tau                 # Stator Length	**GOOD**

ncos1 = np.empty([1,])

R1 = np.empty([0,])
a1 = np.empty([0,])
b1 = np.empty([0,])
X1 = np.empty([0,])
Xm = np.empty([0,])
R2 = np.empty([0,])

Z = np.empty([0,])

I1 = np.empty([0,])
I2 = np.empty([0,])
Im = np.empty([0,])

Fs = np.empty([0,])
diff = np.empty([0,])


for i in range(1,30):
	
	N1 = p*q1*i
	ncos0 = 0.2
	ncos = ncos1[i] = 1
	
	R1 = np.append(i,0)
	a1 = np.append(i,0)
	b1 = np.append(i,0)
	X1 = np.append(i,0)
	Xm = np.append(i,0)
	R2 = np.append(i,0)

	Z = np.append(i,0)

	I1 = np.append(i,0)
	I2 = np.append(i,0)
	Im = np.append(i,0)

	Fs = np.append(i,0)
	diff = np.append(i,0)
	
	while np.abs(ncos0 - ncos1[i]) > 0.0001:
		print 'i:',i
		I1prime = (Fsprime*Vr)/(m*V1*ncos0)
		Aw = I1prime/J1
#		As = (10*i*Aw)/7
		As = (10*i*Aw)/7
		ws = Lambda/2
		wt = ws
		hs = As/ws
		go = gm + d
		gamma = (4/np.pi)*(((ws/(2*go))*np.arctan(ws/(2*go))) - np.log(np.sqrt(1 + ((ws/(2*go))**2))))
		kc = Lambda/(Lambda - gamma*go)
		ge = kc*go
		kw = np.sin(np.pi/(2*m))/(q1*np.sin(np.pi/(2*m*q1)))
		G = 2*mu0*f*tau**2/(np.pi*(rhor/d)*ge)     # Goodness factor for in vacuum
#		G = 2*mu*f*tau**2/(np.pi*(rhor/d)*ge)     # Goodness factor
		a = np.pi/2
		ae = a + ge/2
		Lce = tau
		beta1 = 1
		lamda_s = (hs*(1+3*beta1))/(12*ws)
		lamda_e = (0.3*(3*beta1-1)) 
		lamda_d = 5*(ge/ws)/(5 + 4*(go/ws))

		#Equivalent Circuit Components
		R1[i] = rhow*(4*a + 2*Lce)*J1*N1/I1prime
		a1[i] = lamda_s*(1 + 3/p) + lamda_d
		b1[i] = lamda_e*Lce
		X1[i] = 8*mu0*np.pi*f*((a1[i]*2*a/q1) + b1[i])*N1**2/p # Leakage Reactance
		Xm[i] = (48*mu0*np.pi*f*ae*kw*N1**2*tau)/(np.pi**2*p*ge)  # Magnetizing Reactance
		R2[i] = Xm[i]/G
#		print 'R1[i]',R1[i],'1j',1j,'X1[i]',X1[i],'R2[i]',R2[i],'Xm[i]',Xm[i]
		Z[i] = R1[i] + 1j*X1[i] + ((1j*R2[i]*Xm[i])/Srated)/((R2[i]/Srated) + 1j*Xm[i])
#		print 'V1',V1,'np.abs(Z[i])',np.abs(Z[i])
		I1[i] = V1/np.abs(Z[i])
		I2[i] = 1j*I1[i]*Xm[i]/(R2[i]/Srated + 1j*Xm[i])
		Im[i] = I1[i] - I2[i]

		#Actual TLIM Thrust
#		print 'm:',m, 'I1[i]', I1[i], 'R2[i]',R2[i],'Srated',Srated,'G',G,'Vs',Vs
		Fs[i] = (m*np.abs(I1[i])**2*R2[i])/(((1/(Srated*G)**2) + 1)*Vs*Srated)   # Eqn 3.51
		diff[i] = Fs[i] - Fsprime
		dmin = np.min(np.abs(diff))
		Pout = Fs*Vr	# Eqn 3.4
		Pin = Pout + m*np.abs(I2[i])**2*R2[i] + m*np.abs(I1[i])**2*R1[i]   # Eqn 3.52
		eta = Pout/Pin
		PF = np.cos(np.angle(Z[i]))
		print 'Fs:',Fs
#		print 'Pout:',Pout
#		print 'Pin:',Pin
		print 'eta:',eta
#		print 'PF:',PF
		ncos1[i] = eta*PF
		ncos0 = (ncos0+ncos1[i])/2


dmin = np.min(np.abs(diff))
print dmin

k = 0
while dmin != np.abs(diff[k]):
	k = k + 1


Nc = k         # Number of turns per slot
N1 = p*q1*Nc   # Number of turns per phase
Fs = Fs[k]     # Estimated thrust based on Nc (N)
I1 = I1[k]     # Actual rated stator RMS current (A)

ncos1 = ncos1[k]

A = {	3:	5.8	,
	4:	5.189	,
	5:	4.62	,
	6:	4.1148	,
	7:	3.66	,
	8:	3.2639	,
	9:	2.9057	,
	10:	2.588   }

#gauge = 0

'''
para_wires = np.empty([10,])
slot_width = np.empty([10,])
tooth_width = np.empty([10,])
min_toothwidth = np.empty([10,])
height_slot = np.empty([10,])
Area_wire = np.empty([10,])
Area_slot = np.empty([10,])
Num_c = np.empty([10,])
Num_1 = np.empty([10,])
Sta_I = np.empty([10,])
gap_e = np.empty([10,])
current_den = np.empty([10,])
height_yoke = np.empty([10,])
final_thrust = np.empty([10,])
output = np.empty([10,])
Input = np.empty([10,])
efficiency = np.empty([10,])
difference = np.empty([10,])
diffmin = np.empty([10,])
'''

#while gauge < 8:
for gauge in range(3,10):
	
#	gauge = gauge + 1
#	gauge = 5
	Np = 0
#	r = 0    # Unused variable
	wt = 1
	wtmin = 0
#	g = 0    # Unused variable
#	r = 0    # redundant
	
	while wt - wtmin > 0.0152:
			
#		r = r + 1    # Unused variable
#		g = g + 1    # Unused variable
		Dw = A[gauge]	# Diameter of selected copper wire (mm)
		Np = Np + 1
		ws = (Dw*10**-3*Np) + 2.2*10**-3	# Eqn 4.18
		wt = Lambda - ws	# Eqn 4.19
		Aw = Np*np.pi/4*Dw**2*1e-6
		As = (10*Nc*Aw)/7
		hs = As/ws
		gm = 0.01
		go = gm + d
		gamma = (4/np.pi)*(((ws/(2*go))*np.arctan(ws/(2*go))) - np.log(np.sqrt(1 + ((ws/(2*go))**2))))
		kc = Lambda/(Lambda - gamma*go)
		ge = kc*go
		G = 2*mu0*f*tau**2/(np.pi*(rhor/d)*ge)     # Goodness factor for in vacuum
#		G = 2*mu*f*tau**2/(np.pi*(rhor/d)*ge)     # Goodness factor
		kw = np.sin(np.pi/(2*m))/(q1*np.sin(np.pi/(2*m*q1)))
		a = np.pi/2
		ae = a + ge/2
		Lce = tau
		beta1 = 1
		lamda_s = (hs*(1 + 3*beta1))/(12*ws)
		lamda_e = (0.3*(3*beta1 - 1))
		lamda_d = 5*(ge/ws)/(5 + 4*(go/ws))

		# Equivalent Circuit Components
		R1 = rhow*(4*a + 2*Lce)*J1*N1/I1prime
		a1 = lamda_s*(1 + 3/p)+lamda_d
		b1 = lamda_e*Lce
		X1 = 8*mu0*np.pi*f*((a1*2*a/q1) + b1)*N1**2/p # Leakage Reactance
		Xm = (48*mu0*np.pi*f*ae*kw*N1**2*tau)/(np.pi**2*p*ge)  # Magnetizing Reactance
		R2 = Xm/G
		Z = R1 + 1j*X1 + (R2/Srated*1j*Xm)/(R2/Srated + 1j*Xm)
		I1 = V1/np.abs(Z)
		I2 = 1j*I1*Xm/(R2/Srated + 1j*Xm)
		Im = I1-I2
		wtmin = 2*np.sqrt(2)*m*kw*N1*np.abs(Im)*mu0*Lambda/(np.pi*p*ge*btmax)
		
		
		hy = 4*np.sqrt(2)*m*kw*N1*np.abs(Im)*mu0*Ls/(np.pi*np.pi*p*p*ge*bymax)  # Yoke height of stator core
#		hy = 4*np.sqrt(2)*m*kw*N1*abs(Im)*mu*Ls/(np.pi*np.pi*p*p*ge*bymax)  # Yoke height of stator core
		para_wires[gauge] = Np
		slot_width[gauge] = ws 
		tooth_width[gauge] = wt
		min_toothwidth[gauge] = wtmin
		height_slot[gauge] = hs
		Area_wire[gauge] = Aw
		Area_slot[gauge] = As
		Num_c[gauge] = Nc
		Num_1[gauge] = N1
		Sta_I[gauge] = I1
		gap_e[gauge] = ge
		current_den[gauge] = np.abs(I1)/Aw
		height_yoke[gauge] = 4*np.sqrt(2)*m*kw*N1*(Im)*mu0*Ls/(np.pi*np.pi*p*p*ge*bymax)
#		height_yoke[gauge] = 4*np.sqrt(2)*m*kw*N1*(Im)*mu*Ls/(np.pi*np.pi*p*p*ge*bymax)
		final_thrust[gauge] = (m*np.abs(I1)**2*R2)/(((1/(Srated*G)**2) + 1)*Vs*Srated)
		output[gauge] = final_thrust[gauge]*Vr
		Input[gauge] = output[gauge] + m*np.abs(I2)**2*R2 + m*np.abs(I1)**2*R1
		efficiency[gauge] = output[gauge]/Input[gauge]
		difference[gauge] = final_thrust[gauge] - Fsprime
		diffmin[gauge] = np.min(np.abs(difference))
'''
		np.insert(para_wires, gauge, Np )
		np.insert(slot_width, gauge, ws )
		np.insert(tooth_width, gauge, wt )
		np.insert(min_toothwidth, gauge, wtmin )
		np.insert(height_slot, gauge, hs )
		np.insert(Area_wire, gauge, Aw )
		np.insert(Area_slot, gauge, As )
		np.insert(Num_c, gauge, Nc )
		np.insert(Num_1, gauge, N1 )
		np.insert(Sta_I, gauge, I1 )
		np.insert(gap_e, gauge, ge )
		np.insert(current_den, gauge, np.abs(I1)/Aw )
		np.insert(height_yoke, gauge, 4*np.sqrt(2)*m*kw*N1*(Im)*mu0*Ls/(np.pi*np.pi*p*p*ge*bymax) )
#		height_yoke[gauge] = 4*np.sqrt(2)*m*kw*N1*(Im)*mu*Ls/(np.pi*np.pi*p*p*ge*bymax)
		np.insert(final_thrust, gauge,  (m*np.abs(I1)**2*R2)/(((1/(Srated*G)**2) + 1)*Vs*Srated) )
		np.insert(output, gauge,  final_thrust[gauge]*Vr )
		np.insert(Input, gauge,  output[gauge] + m*np.abs(I2)**2*R2 + m*np.abs(I1)**2*R1 )
		np.insert(efficiency, gauge, output[gauge]/input[gauge] )
		np.insert(difference, gauge, final_thrust[gauge] - Fsprime )
		np.insert(diffmin, gauge, np.min(np.abs(difference)) )
'''

kk = np.min(diffmin)
jj = 1

while kk != np.abs(diffmin(jj)):
	jj = jj + 1



best_wiregauge = A(jj,1)

#$$$ To Generate the Characteristic curves $$$

# vel_sta = 17.22    # Only applies for Vr = 15.5 & Srated = 10#
# slip = 0.1         # Only applies for Srated = 10#
vel_sta = Vs
# slip = Srated    # commented out, b/c variable gets overwritten anyway
e = 1

for slip in range(0.000001,1,0.01):
    
    vel_rot[e] = vel_sta*(1 - slip)
    if np.abs(Vr - vel_rot[e])/Vr < 0.01:
        n_Vr = e   # Index for where v = Vr
    end
    impz[e] = R1 + 1j*X1 + (R2/slip*1j*Xm)/(R2/slip + 1j*Xm)
    i1[e] = V1/np.abs(impz[e])
    i2[e] = 1j*i1[e]*Xm/(R2/slip + 1j*Xm)
    im[e] = i1[e] - i2[e]
    Force[e] = (m*(np.abs(i1[e]))**2*R2)/(((1/(slip*G)**2) + 1)*vel_sta*slip)
    out_pow[e] = Force[e]*vel_rot[e]
    in_pow[e] = out_pow[e] + m*np.abs(i2[e])**2*R2 + m*np.abs(i1[e])**2*R1
    eff[e] = out_pow[e]/in_pow[e]
    e = e + 1
    
end

## Miscellaneous calcs (not included in original code)
# Physical properties
rhoiron = 7870                     # Density of iron (kg/m**3)
rhocopper = 8960                   # Density of copper (kg/m**3)

# Assumptions & Dummy values
lce = 0.1144                       # Length of end connection

# Ammount of material required for construction of one SLIM stator
lw = 2*(Ws + lce)*N1               # Length of one turn of copper winding inside a stator slot
Tlw = m*Np*lw                      # Length of copper wire required for stator windings

Vyoke = Ls*Ws*hy                   # Volume of iron required for stator yoke
Vtooth = Ws*wt*hs                  # Volume of iron required for stator tooth
Vteeth = m*p*q1*Ws*wt*hs           # Volume of iron required for stator teeth
Viron = Ws*(Ls*hy + m*p*q1*wt*hs)  # Total volume of iron required
Vcopper = Tlw*np.pi*(Dw*10**-3/2)**2	# Volume of copper used (m**3)

Wiron = rhoiron*Viron              # Total weight of iron required
Wcopper = rhocopper*Vcopper        # Total weight of copper required
Wstator = Wiron + Wcopper          # Total weight of stator


## Generate table of outputs 

table = [
	#VariableNames 
	[
            'Rated Slip'
            'Yoke density'
            'Tooth density'
            'Core Width'
            'SLIM Synchronous velocity'
            'Rotor velocity'
            'No. of Poles'
            'Pole np.pitch'
            'Slot np.pitch'
            'Stator length'
            '"Target" Thrust'
            'No. of turns per slot'
            'No. of turns per phase'
            'Copper wire size in winding'
            'Diameter of selected copper wire (mm)'
            'Parallel wires'
            'Slot width'
            'Tooth width'
            'Minimum tooth width'
            'Slot depth'
            'Stator core yoke height'
            'Actual thrust at specified Vr'
            'Output power at specified Vr'
            'Input power at specified Vr'
            'Stator efficiency at specified Vr'
            'Actual rated stator RMS current'
            'Actual stator current density'
            'Total length of copper wire'
            'Total weight of copper wire'
            'Iron core weight'
            'Total weight of one stator unit'
						],

	# Variable
	[   
	    'Srated'
	    'bymax'
	    'btmax'
	    'Ws'
	    'Vs'
	    'Vr'
	    'p'
	    'tau'
	    'lamda'
	    'Ls'
	    'Fsprime'
	    'Nc'
	    'N1'
	    'gauge'
	    'Dw'
	    'Np'
	    'ws'
	    'wt'
	    'wtmin'
	    'hs'
	    'hy'
	    'Fs(Vr)'
	    'Pout(Vr)'
	    'np.pin(Vr)'
	    'eta(Vr)'
	    'I1(Vr)'
	    'J1'
	    'Tlw'
	    'Wcopper'
	    'Wiron'
	    'Wstator'
			],

	#Value	
	[	
		Srated,
		bymax,
		btmax,
		Ws,
		Vs,
		Vr,
		p,
		tau,
		Lambda,
		Ls,
		Fsprime,
		Nc,
		N1,
		gauge,
		Dw,
		Np,
		ws,
		wt,
		wtmin,
		hs,
		hy,
		Force(Vr),
		out_pow(Vr),
		in_pow(Vr),
		eff(Vr),
		I1,
		J1,
		Tlw,
		Wcopper,
		Wiron,
		Wstator
		    	],

	#Units
	[   
		'-'
		'Tesla'
		'Tesla'
		'm'
		'm/s'
		'm/s'
		'-'
		'm'
		'm'
		'm'
		'N'
		'-'
		'-'
		'AWG'
		'mm'
		'-'
		'm'
		'm'
		'm'
		'm'
		'm'
		'N'
		'W'
		'W'
		'-'
		'A'
		'A/m**2'
		'm'
		'kg'
		'kg'
		'kg'
			],

	#Dependencies
	[   
		'Design Parameter/Constant'
		'Design Parameter/Constant'
		'Design Parameter/Constant'
		'Design Parameter/Constant'
		'Vr,Srated'
		'Design Parameter/Constant'
		'Design Parameter/Constant'
		'Vs,f'
		'tau,m,q1'
		'p,tau'
		'Design Parameter/Constant'
		'k'
		'p,q1,Nc'
		'(see loop counter)'
		'gauge'
		'(see loop counter)'
		'Dw,Np'
		'Lambda,ws'
		'm,kw,N1,Im,mu0,Lambda,p,ge,btmax'
		'Asws'
		'm,kw,N1,Im,mu0,Ls,np.pi,p,ge,bymax'
		'm,i1[e],R2,slip,G,vel_sta'
		'Force[e],vel_rot[e]'
		'out_pow[e],m,i2[e],R2,m,i1[e],R1'
		'out_pow[e],in_pow[e]'
		'V1,Z[i]'
		'Design Parameter/Constant'
		'm,Np,lw'
		'rhocopper,Vcopper'
		'rhoiron,Viron'
		'Wcopper,Wiron'
        		],
			]


import csv

# write it
with open('SLIM_results.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow(r) for r in table]

# read it
with open('SLIM_results.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    table = [[int[e] for e in r] for r in reader]

'''
T = table(Variable,Value,Units,Dependencies,line_reference,'RowNames',VariableNames)

# Write table of output parameters to csv
formatSpec = 'SLIM_case_no_#0.f.csv'
filename = sprintf(formatSpec,designno)
writetable(T,filename)
'''

'''

## Graph Thrust and Efficiency

figure(1)
plt.plot(vel_rot,Force,[Vr, Vr],[0,Fs])
#hold on
#grid on
#grid minor
ylabel('Target Thrust, Fs (N)')
xlabel('Rotor Velocity, Vr (m/s)')
# plt.plot([15.5 15.5],[0,Fs])  # Only applies for Vr = 15.5
#plt.plot([Vr, Vr],[0,Fs])
#hold on
# plt.plot([0 15.5],[Fs Fs])  # Only applies for Vr = 15.5
plt.plot([0 Vr],[Fs Fs])
#hold on
title(['Force vs. Velocity for design case no. ' num2str(designno)])
legend('Actual Force','Target Velocity','Target Force')

figure(2)
plt.plot(vel_rot,eff*100,'green')
#hold on
# plt.plot([15.5 15.5],[0 eta*100])  # Only applies for Vr = 15.5
plt.plot([Vr Vr],[0 eta*100])
#hold on
# plt.plot([0 Vr],[eta*100,eta*100])  # Only applies for Vr = 15.5
plt.plot([0 Vr],[eta*100,eta*100])
hold on
grid on
grid minor
ylabel('Efficiency (#)')
xlabel('Rotor Velocity, Vr (m/s)')
title(['Efficiency vs. Velocity for design case no. ' num2str(designno)])
legend('Actual Efficiency','Target Velocity','Ideal Efficiency')

## Write table of rotor velocity dependent outputs to csv
VariableNames = {   
                    'Rotor_Velocity',
                    'Impedance',
                    'I1',
                    'I2',
                    'Im',
                    'Thrust',
                    'Output_power',
                    'Input_power',
                    'Efficiency'
                }

size(VariableNames)

T = table(vel_rot', impz', i1', i2', im', Force', out_pow', in_pow', eff','VariableNames',VariableNames)

formatSpec = 'LIM rotor vel dependent outputs, case no. #0.f.csv'
filename = sprintf(formatSpec,designno)
writetable(T,filename)

## Assign Operating parameters
# Reduced voltage and frequency with reduced rated rotor velocity
# run simulation at design parameters without reloonp.ping to find 
# new values except calculate at the three changed constants

operatingno = 13   # Set Design Case no.

switch operatingno
    
    case 10 # 
		factor = .001
		Vline = Vline * factor
		f = f * factor
		Vr = Vr * factor

	case 11 # 
		factor = .1
		Vline = Vline * factor
		f = f * factor
		Vr = Vr * factor
		
	case 12 #
		factor = .2
		Vline = Vline * factor
		f = f * factor
		Vr = Vr * factor
		
	case 13 #
		factor = .3
		Vline = Vline * factor
		f = f * factor
		Vr = Vr * factor
		
	case 14 #
		factor = .4
		Vline = Vline * factor
		f = f * factor
		Vr = Vr * factor
		
	case 19 #
		factor = .9
		Vline = Vline * factor
		f = f * factor
		Vr = Vr * factor
        
end
		
## Calculate output performance per design case and operating conditions

# Data from the PCP design procedure
V1 = Vline/np.sqrt(3)         # Rated primary RMS - Eqn 4.16
Vs = Vr/(1 - Srated)       # Sychronous velocity (m/s)	**GOOD**
tau = Vs/(2*f)             # Pole np.pitch	**GOOD**
Lambda = tau/(m*q1)        # Slot np.pitch    **GOOD**
Ls = p*tau                 # Stator Length	**GOOD**

Sta_I[gauge] = I1
current_den[gauge] = abs(I1)/Aw

vel_sta = Vs
# slip = Srated    # commented out, b/c variable gets overwritten anyway
e = 1

for slip = 0.000001:0.01:1
    
    vel_rot[e] = vel_sta*(1 - slip)
    impz[e] = R1 + 1j*X1 + (R2/slip*1j*Xm)/(R2/slip + 1j*Xm)
    i1[e] = V1/abs(impz[e])
    i2[e] = 1j*i1[e]*Xm/(R2/slip + 1j*Xm)
    im[e] = i1[e] - i2[e]
    Force[e] = (m*(abs(i1[e]))**2*R2)/(((1/(slip*G)**2) + 1)*vel_sta*slip)
    out_pow[e] = Force[e]*vel_rot[e]
    in_pow[e] = out_pow[e] + m*abs(i2[e])**2*R2 + m*abs(i1[e])**2*R1
    eff[e] = out_pow[e]/in_pow[e]
    e = e + 1
    
end

## Graph Thrust and Efficiency for operating conditions

figure(3)
plt.plot(vel_rot,Force, [Vr, Vr],[0,Fs], [0 Vr],[Fs, Fs])
#plot(vel_rot,Force,'green')
#hold on
#grid on
#grid minor
plt.ylabel('Target Thrust, Fs (N)')
plt.xlabel('Rotor Velocity, Vr (m/s)')
# plt.plot([15.5 15.5],[0,Fs])  # Only applies for Vr = 15.5
#plt.plot([Vr Vr],[0,Fs])
#hold on
# plt.plot([0 15.5],[Fs Fs])  # Only applies for Vr = 15.5
#plt.plot([0 Vr],[Fs Fs])
#hold on
plt.title(['Force vs. Velocity Outputs for design case no. ' num2str(designno) ' and operating case no. ',num2str(operatingno)])
plt.legend('Actual Force','Target Velocity','Target Force')

figure(4)
#plt.plot(vel_rot,eff*100,'green')
plt.plot(vel_rot,eff*100, [Vr Vr],[0 eta*100], [0 Vr],[eta*100,eta*100])
#hold on
# plt.plot([15.5 15.5],[0 eta*100])  # Only applies for Vr = 15.5
#plt.plot([Vr Vr],[0 eta*100])
#hold on
# plt.plot([0 Vr],[eta*100,eta*100])  # Only applies for Vr = 15.5
#plt.plot([0 Vr],[eta*100,eta*100])
#hold on
#grid on
#grid minor
plt.ylabel('Efficiency (#)')
plt.xlabel('Rotor Velocity, Vr (m/s)')
plt.title(['Efficiency vs. Velocity Outputs for design case no. ' num2str(designno) ' and operating case no. ',num2str(operatingno)])
plt.legend('Actual Efficiency','Target Velocity','Ideal Efficiency')

plt.show()
'''

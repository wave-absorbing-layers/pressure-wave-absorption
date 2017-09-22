# -*- coding: utf-8 -*-
import sys
import math
import cmath
#import pylab

try:
    import pylab
    pylab_available = "YES"
except:
    pylab_available = "NO"


####### PROGRAM DESCRIPTION ###########################
#
### Theory_1D.py ###
# written by Robinson Peric 
# at Hamburg University of Technology (TUHH)
# version: 22 september 2017
#
### What the code does ###
# The code computes analytical predictions for the reflection coefficient
# for 1D-pressure wave propagation when using 
# absorbing layer approaches as described in [1]
#
### Please note ###
# In [1] it is found that for many  2D- and 3D-flow simulations, 1D-theory predictions (as this code provides) are sufficiently accurate. 
# Recommendations are given in [1] how to tune the absorbing layer parameters depending on the waves. 
# Although 1D-theory predictions are often quite accurate [1], please keep in mind that every theory has its limitations! Tuning the absorbing layer parameters according to the present theory does NOT guarantee that the actual reflection coefficient in the simulation will equal the prediction. Other mechanisms can also lead to undesired wave reflections. Further, the behavior of absorbing layers for highly non-linear waves, such as distorted or shock waves, is not fully understood so far. See [1] for a detailed discussion.
# See also the manual distributed along with the code.
#
### References ###
# [1] Robinson Perić. Analytical Prediction of Reflection Coefficients for Wave Absorbing Layers. Preprint, arXiv:1705.06937 [physics.flu-dyn], 2017.
#
### How to use the program ###
# Call program like this:
# python ./Theory_1D.py
#
### Requirements ###
# The programming language python version 3.0 or higher (https://www.python.org/downloads/) must be installed. 
# It is recommended to also have matplotlib (https://matplotlib.org/users/installing.html) installed. Then the code will plot results beautifully.
#
# Please report errors to robinson.peric@tuhh.de
# 
#####################################################




####### ENTER USER ARGUMENTS ###########################

foo = input("\n\n\nPlease enter wave period (s):\n")
T = float(foo)
foo = input("\nPlease enter speed of sound (m/s):\n")
L = float(foo) * T 	# convert speed of sound to wavelength
foo = input("\nPlease enter layer thickness per wavelength (enter 2 if thickness = 2 * wavelength):\n")
layerThickness = L*float(foo)
print("\nAvailable blending functions are:")
print("Constant blending: 1")
print("Linear blending: 2")
print("Quadratic blending: 3")
print("Cosine**2 blending: 4")
print("Exponential blending: 5")
print("Custom blending: 6")
foo = input("\nPlease enter number of desired blending function:\n")
blend = float(foo)


### default settings
csv_file_separator = " "	# alt. "," or ";"
factor =1.05				# resolution of plot: very fine (1.01), very coarse (4.0)
gammaMin=10**-4/T			# range of forcing strength gamma for which reflection coefficient is computed
gammaMax=10**7/T			# -"-
dampres=200					# number of piece-wise constant blending zones, into which absorbing layer is subdivided
Lx = 2.0*L 	# for plots: domain size outside the forcing layer







###### VARIABLE DECLARATIONS AND INITIALIZATIONS ###########################

# calculate wave parameters 
w = 2 * math.pi / T
wavenumber_k = 2 * math.pi / L
c = L / T	# phase velocity

# initialize vectors
xd = [0.0] * ( int(dampres) + 2 )		 	# initialize vector for thicknesses of layers, xd[1] is Lx
k = [0.0] * ( int(dampres) + 2 )		 	# initialize vector for wave number, k[1] is 2*pi/L
Ct = [0.0] * ( int(dampres) + 2 )		 	# initialize vector for transmission coefficient, Ct[1] is coefficient for entrance to dampign zone
beta = [0.0] * ( int(dampres) + 2 )			# initialize vector needed to compute Cr
Cr = [0.0] * ( int(dampres) + 2 )		 	# initialize vector for reflection coefficient, Cr[1] is the global reflection coefficient (i.e. for the entrance to the forcing layer)

# write initial conditions
dx = layerThickness / dampres	# thickness of single 'cell' within damping layer
xd[1] = Lx
for i in range(2,dampres+2,1):
	xd[i] = dx
k[0] = wavenumber_k
k[1] = wavenumber_k
Ct[0] = 1.0
Cr[0] = 0.0
Ct[dampres+1] = 0.0
Cr[dampres+1] = 1.0




###### FUNCTION DECLARATIONS ###########################

def abs(x):
	return math.sqrt( (x.real)**2 + (x.imag)**2 )	# return absolute value of complex number x

# evaluate blending function
def b(x):
	if (x >= Lx):

		# constant	
		if (blend == 1):
			return 1.0	# constant blending

		# linear
		elif (blend == 2):
			return  ( x - Lx ) / ( ( Lx + layerThickness ) - Lx )

		# quadratic
		elif (blend == 3):
			return  ( ( x - Lx ) / ( ( Lx + layerThickness ) - Lx ) )**2

		# cosine-squared
		elif (blend == 4):
			return math.cos( 0.5 * math.pi + 0.5 * math.pi * ( x - Lx ) / ( ( Lx + layerThickness ) - Lx ) )**2

		# exponential
		elif (blend == 5):
			return (math.exp( ( ( x - Lx ) / ( ( Lx + layerThickness ) - Lx )  )**2 ) -1) / ( math.exp(1) - 1 )

		# custom
		elif (blend == 6):
			return 42 # replace "42" by a user-defined custom blending function, which lies in [0,1]
		
		else:
			print("ERROR: Specified blending function not found.")
	
	else:
		return 0.0

def fbeta(i):
    return ( ( 1.0 + Cr[i] * cmath.exp( 1j * k[i] * xd[i] * 2.0 )  ) / ( 1.0 -  Cr[i] * cmath.exp( 1j * k[i] * xd[i] * 2.0 ) ) )

def fCt(i):
    return ( 1.0 - Cr[i] ) / (  1.0 - Cr[i+1] * cmath.exp( 1j * k[i+1] * xd[i+1] * 2.0 )  ) 

def fCr(i):
    return ( -k[i] + k[i+1] * beta[i+1] ) / ( k[i] + k[i+1] * beta[i+1] )






###### MAIN ##############################################
gamma = gammaMin
G=[]
CR=[]

# create file
f = open("C_R.csv", 'w')
f.close()

# calculate C_R for different gamma values
while (gamma < gammaMax):
	tmp = Lx
	for i in range(2,dampres+2,1):
		k[i] = cmath.sqrt(  (w**2)/ ( c**2 ) + (1j * w * gamma * b( tmp + 0.5 * xd[i] ) ) / ( (c)**2 ))
		tmp += xd[i]

	# calculate transmission and reflection coefficients
	for i in range(dampres,0,-1):   #go cell-by-cell from end of the layer to layer entrance. E.g. for 3 zones, start with zone 3, then 2, then 1.
		beta[i+1] = fbeta(i+1)
		Cr[i] = fCr(i)
		Ct[i] = fCt(i)

	# write reflection coefficients to file
	f = open("C_R.csv", 'a')
	f.write(str(gamma) + csv_file_separator + str(abs(Cr[1])) +  "\n")
	G.append(gamma)
	CR.append(abs(Cr[1]))
	f.close()

	gamma = gamma * factor


# if pylab and matplotlib are installed, open window and plot C_R(gamma)
if (pylab_available == "YES"):
	pylab.figure(figsize=(20,8))
	pylab.xlabel('$\gamma\ (1/\mathrm{s})$ ', {'fontsize': 20})
	pylab.ylabel('$C_R$', {'fontsize': 20})
	pylab.semilogx()
	pylab.semilogy()
	pylab.plot(G,CR, linewidth=1.5)
	pylab.show()

print("\n\nProgram finished. \nReflection coefficients were written to file C_R.csv\n")

sys.exit(0)


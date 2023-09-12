#----------------------------------------------------------------------------------------
#  ______ ____   _____   _  _     _______     _________ _    _  ____  _   _ 
# |  ____/ __ \ / ____| | || |   |  __ \ \   / /__   __| |  | |/ __ \| \ | |
# | |__ | |  | | (___   | || |_  | |__) \ \_/ /   | |  | |__| | |  | |  \| |
# |  __|| |  | |\___ \  |__   _| |  ___/ \   /    | |  |  __  | |  | | . ` |
# | |___| |__| |____) |    | |   | |      | |     | |  | |  | | |__| | |\  |
# |______\____/|_____/     |_|   |_|      |_|     |_|  |_|  |_|\____/|_| \_|
#
#----------------------------------------------------------------------------------------
#
#
#   Technische Universiteit Delft - TUDelft (2022)
#
#   Master of Science in Chemical Engineering
#
#   This code was developed and tested by Elia Ferretti
#
#   You can redistribute the code and/or modify it
#   Whenever the code is used to produce any publication or document,
#   reference to this work and author should be reported
#   No warranty of fitness for a particular purpose is offered
#   The user must assume the entire risk of using this code
#
#
#----------------------------------------------------------------------------------------     

# import the library
from thermodynamicEOS import EOS

#----------------------------------------------------------------------------------------  
# AIR MIXTURE
#----------------------------------------------------------------------------------------  

# define the critical constants for the species in the mixture
#   - [0] oxygen O2
#   - [1] nitrogen N2
Tc = [154.5,    126.19]         # [K]
pc = [50.43e5,  33.98e5]        # [Pa]
w =  [0.025,    0.037]          # [-]

# define conditions and mixture property
T = 300                         # [K]   - (float) absolute temperature
p = 10e5                        # [Pa]  - (float) absolute pressure
x = [0.21, 0.79]                # [-]   - (list)  molar fractions composition
state = "V"                     # [-]   - (char)  mixture state (either "V" or "L")

# testing the checks
p,T,x,state,index,fatalError = EOS.check("C",0,0,pc,Tc,[0.2,0.7],w)

print("\n")
print("hR_PRmix =\t",EOS.hR_PRmix(Tc, pc, w, T, p, x, state),"\t[J/mol]")
print("hR_RKSmix =\t",EOS.hR_RKSmix(Tc, pc, w, T, p, x, state),"\t[J/mol]")
print("hR_RKmix =\t",EOS.hR_RKmix(Tc, pc, T, p, x, state),"\t[J/mol]")
print("hR_VdWmix =\t",EOS.hR_VdWmix(Tc, pc, T, p, x, state),"\t[J/mol]")


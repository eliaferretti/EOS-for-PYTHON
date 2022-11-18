Extended library containing functions for calculating thermodynamic quantities with various cubic EoS + Virial

EoS prerequisites Declare global variables:
- Tc = Vector of critical temperatures [K]
- pc = Vector of critical pressures [Pa]
- w = Vector of Pitzer acentric factors [-]

Meaning of input terms:
- T/Temp = Temperature [K]
- p/press = Pressure [Pa]
- x = phase composition (vector of molar fractions)
- state = physical state ("L" or "V")
- i/index = index of the compound

Notes on output terms:
- zeta [-]
- phi [-]
- residual enthalpy [J/mol]
- gamma [-]

Notes - The functions use the units of measurement of the International System - All functions have been tested on the numerical results of exercises proposed in the text: "Fondamenti di Termodinamica dell’Ingegneria Chimica" - R.Rota

For any problems or bugs, please contact the developer at: eliaferretti@outlook.it

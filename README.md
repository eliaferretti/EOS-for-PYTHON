## EOS for PYTHON ##

Extended library containing functions for calculating thermodynamic quantities with various cubic **Equations of State**.

### Meaning of input terms: ###
- `Tc` = vector of critical temperatures $\left[K\right]$
- `pc` = vector of critical pressures $\left[Pa\right]$
- `w` = vector of Pitzer acentric factors $\left[-\right]$
- `T`/`Temp` = temperature $\left[K\right]$
- `p`/`press` = pressure $\left[Pa\right]$
- `x` = mixture composition, vector of molar fractions
- `state` = physical state ("L" or "V")
- `i`/`index` = index of the compound

### Notes on output terms: ###

- `zeta` = compressibility factor $\left[-\right]$
- `phi` = fugacity coefficient $\left[-\right]$
- `residual enthalpy` $\left[\frac{J}{mol}\right]$

### Additional notes: ###
- The functions use the units of measurement of the International System
- All functions have been tested on the numerical results of exercises proposed in the text: "*Fondamenti di Termodinamica dell’Ingegneria*" - R.Rota
- For any problems or bugs, please contact the developer at: **eliaferretti@outlook.it**

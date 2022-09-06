#
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
#---------------------------------------------------------------------------------------                                                                         
                                                                           
def zeta_VdW(Tc,pc,T,pressure,state):
    import math
    R = 8.3144621
    
    #Tr = T / Tc
    a_min = 0.421875 * (R * Tc) ** (2) / pc
    b_min = 0.125 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)
    pig = 3.14159265353589
    
    beta = - 1 - b
    gamma = a
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta_VdW = aus
    return zeta_VdW
    
def zeta_RK(Tc,pc,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    k = 1/(Tr) ** (0.5)
    a_min = 0.42748 * (R * Tc) ** (2) * k / pc
    b_min = 0.08664 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)
    pig = 3.14159265353589
    
    beta = - 1
    gamma =  a - b - (b) ** (2)
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta_RK = aus
    return zeta_RK

def zeta_RKS(Tc,pc,w,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    S = 0.48 + 1.574 * w - 0.176 * (w) ** (2)
    k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
    a_min = 0.42748 * (R * Tc) ** (2) * k / pc
    b_min = 0.08664 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)
    pig = 3.14159265353589
    
    beta = - 1
    gamma = a - b - (b) ** (2)
    delta = - a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta_RKS = aus
    return zeta_RKS
    
def zeta_PR(Tc,pc,w,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    S = 0.37464 + 1.54226 * w - 0.26992 * w ** (2)
    k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
    a_min = 0.45724 * (R * Tc) ** (2) * k / pc
    b_min = 0.0778 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)
    pig = 3.14159265353589
    
    beta = - 1 + b
    gamma = a - 2 * b - 3 * (b) ** (2)
    delta = -a * b + (b) ** (2) + (b) ** (3)
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta_PR = aus
    return zeta_PR


   
def phi_VdW(Tc,pc,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    a_min = 0.421875 * (R * Tc) ** (2) / pc
    b_min = 0.125 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)
    pig = 3.14159265353589
    
    beta = - 1 - b
    gamma = a
    delta = - a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta = aus
    aus = zeta - 1 - a / zeta - math.log( zeta - b )
    phi_VdW = math.exp(aus)
    return phi_VdW
    
def phi_RK(Tc,pc,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    k = 1/(Tr) ** (0.5)
    a_min = 0.42748 * (R * Tc) ** (2) * k / pc
    b_min = 0.08664 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)
    pig = 3.14159265353589
    
    beta = - 1
    gamma = a - b - (b) ** (2)
    delta = - a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta = aus
    aus = zeta - 1 - a * math.log((zeta + b) / zeta) / b - math.log(zeta - b)
    phi_RK = math.exp(aus)
    return phi_RK

def phi_RKS(Tc,pc,w,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    S = 0.48 + 1.574 * w - 0.176 * (w) ** (2)
    k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
    a_min = 0.42748 * (R * Tc) ** (2) * k / pc
    b_min = 0.08664 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)
    pig = 3.14159265353589
    
    beta = - 1
    gamma = a - b - (b) ** (2)
    delta = - a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta = aus
    aus = zeta - 1 - a * math.log((zeta + b) / zeta) / b - math.log(zeta - b)
    phi_RKS = math.exp(aus)
    return phi_RKS
    
def phi_PR(Tc,pc,w,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    S = 0.37464 + 1.54226 * w - 0.26992 * w ** (2)
    k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
    a_min = 0.45724 * (R * Tc) ** (2) * k / pc
    b_min = 0.0778 * R * Tc / pc

    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)
    pig = 3.14159265353589
    
    beta = - 1 + b
    gamma = a - 2 * b - 3 * (b) ** (2)
    delta = -a * b + (b) ** (2) + (b) ** (3)
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta = aus
    aus = zeta - 1 - a * math.log((zeta + b) / zeta) / b - math.log(zeta - b)
    phi_PR = math.exp(aus)
    return phi_PR
    


def hR_VdW(Tc,pc,T,pressure,state):
    import math
    
    R = 8.3144621
    
    #Tr = T / Tc
    a_min = 0.421875 * (R * Tc) ** (2) / pc
    b_min = 0.125 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)

    pig = 3.14159265353589
    
    beta = - 1 - b
    gamma = a
    delta = - a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta = aus
    
    hR_VdW = (zeta - 1 - a / zeta)*R*T
    return hR_VdW

def hR_RK(Tc,pc,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    k = 1/(Tr) ** (0.5)
    a_min = 0.42748 * (R * Tc) ** (2) * k / pc
    b_min = 0.08664 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)

    pig = 3.14159265353589
    
    beta = - 1
    gamma = a - b - (b) ** (2)
    delta = - a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta = aus
    
    hR_RK = (zeta - 1 - 1.5 * a / b * math.log((zeta + b)/zeta))*R*T
    return hR_RK

def hR_RKS(Tc,pc,w,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    S = 0.48 + 1.574 * w - 0.176 * (w) ** (2)
    k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
    a_min = 0.42748 * (R * Tc) ** (2) * k / pc
    b_min = 0.08664 * R * Tc / pc
    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)

    pig = 3.14159265353589
    
    beta = - 1
    gamma = a - b - (b) ** (2)
    delta = - a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta = aus
    
    e = S * (Tr / k)**(0.5)
    hR_RKS = (zeta - 1 - ( 1 + e ) * a / b * math.log((zeta + b)/zeta))*R*T
    return hR_RKS 
    
def hR_PR(Tc,pc,w,T,pressure,state):
    import math
    
    R = 8.3144621
    
    Tr = T / Tc
    
    S = 0.37464 + 1.54226 * w - 0.26992 * w ** (2)
    k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
    a_min = 0.45724 * (R * Tc) ** (2) * k / pc
    b_min = 0.0778 * R * Tc / pc

    a = a_min * pressure / (R * T) ** (2)
    b = b_min * pressure / (R * T)

    pig = 3.14159265353589
    
    beta = - 1 + b
    gamma = a - 2 * b - 3 * (b) ** (2)
    delta = - a * b + (b) ** (2) + (b) ** (3)
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)
        coeff1 = 1
        coeff2 = 1
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1
        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3
    zeta = aus
    
    e = S * (Tr / k)**(0.5)
    hR_PR = (zeta - 1 - ( 1 + e ) * a / b / (2)**(1.5) * math.log((zeta + b * (1 + (2)**(0.5)))/(zeta + b * (1 - (2)**(0.5)))))*R*T
    return hR_PR
    
  
  
def zeta_VdWmix(Tc,pc,T,pressure,x,state):
    
    import math
    
    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        #Tr = T / Tc[i]
        Ac[i] = 0.421875 * (R * Tc[i]) ** (2) / pc[i]
        Bc[i] = 0.125 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)

    pig = 3.14159265353589
	
    beta = - 1 - b
    gamma = a
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus

    zeta_VdWmix = zeta
    return zeta_VdWmix

def zeta_RKmix(Tc,pc,T,pressure,x,state):
    
    import math
    
    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        Tr = T / Tc[i]
        k = 1 / (Tr)**(0.5)
        Ac[i] = 0.42748 * (R * Tc[i]) ** (2) * k / pc[i]
        Bc[i] = 0.08664 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)

    pig = 3.14159265353589
	
    beta = -1
    gamma = a - b - (b) ** (2)
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus

    zeta_RKmix = zeta
    return zeta_RKmix
    
def zeta_RKSmix(Tc,pc,w,T,pressure,x,state):
    
    import math
    
    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        S = 0.48 + 1.574 * w[i] - 0.176 * (w[i]) ** (2)
        Tr = T / Tc[i]
        k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
        Ac[i] = 0.42748 * (R * Tc[i]) ** (2) * k / pc[i]
        Bc[i] = 0.08664 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)

    pig = 3.14159265353589
	
    beta = -1
    gamma = a - b - (b) ** (2)
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus

    zeta_RKSmix = zeta
    return zeta_RKSmix

def zeta_PRmix(Tc,pc,w,T,pressure,x,state):
    
    import math
    
    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        S = 0.37464 + 1.54226 * w[i] - 0.26992 * (w[i]) ** (2)
        Tr = T / Tc[i]
        k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
        Ac[i] = 0.45724 * (R * Tc[i]) ** (2) * k / pc[i]
        Bc[i] = 0.0778 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)

    pig = 3.14159265353589
	
    beta = - 1 + b
    gamma = a - 2*b - 3*(b) ** (2)
    delta = -a * b + (b)**(2) + (b)**(3)
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus

    zeta_PRmix = zeta
    return zeta_PRmix



def phi_VdWmix(Tc,pc,index,T,pressure,x,state):
    
    import math
    index = index - 1
    
    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        #Tr = T / Tc[i]
        Ac[i] = 0.421875 * (R * Tc[i]) ** (2) / pc[i]
        Bc[i] = 0.125 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)
    Ai = Ac[index] * pressure / (R * T) ** 2
    Bi = Bc[index] * pressure / (R * T)

    pig = 3.14159265353589
	
    beta = - 1 - b
    gamma = a
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus
    
    aus = Bi/(zeta-b)-2*(Ai*a)**(0.5)/zeta-math.log(zeta-b) 
    
    phi_VdWmix = math.exp(aus)

    return phi_VdWmix

def phi_RKmix(Tc,pc,index,T,pressure,x,state):
    
    import math
    index = index - 1
    
    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        Tr = T / Tc[i]
        k = 1 / (Tr)**(0.5)
        Ac[i] = 0.42748 * (R * Tc[i]) ** (2) * k / pc[i]
        Bc[i] = 0.08664 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)
    Ai = Ac[index] * pressure / (R * T) ** 2
    Bi = Bc[index] * pressure / (R * T)

    pig = 3.14159265353589
	
    beta = -1
    gamma = a - b - (b) ** (2)
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus
    
    aus = Bi * (zeta - 1) / b + a * (Bi / b - 2 * (Ai / a) ** (0.5)) * math.log((zeta + b) / zeta) / b - math.log(zeta - b) 
    
    phi_RKmix = math.exp(aus)

    return phi_RKmix
    
def phi_RKSmix(Tc,pc,w,index,T,pressure,x,state):
    
    import math
    index = index - 1
    
    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        S = 0.48 + 1.574 * w[i] - 0.176 * (w[i]) ** (2)
        Tr = T / Tc[i]
        k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
        Ac[i] = 0.42748 * (R * Tc[i]) ** (2) * k / pc[i]
        Bc[i] = 0.08664 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)
    Ai = Ac[index] * pressure / (R * T) ** 2
    Bi = Bc[index] * pressure / (R * T)

    pig = 3.14159265353589
	
    beta = -1
    gamma = a - b - (b) ** (2)
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus
    
    aus = Bi * (zeta - 1) / b + a * (Bi / b - 2 * (Ai / a) ** (0.5)) * math.log((zeta + b) / zeta) / b - math.log(zeta - b) 
    
    phi_RKSmix = math.exp(aus)

    return phi_RKSmix
    
def phi_PRmix(Tc,pc,w,index,T,pressure,x,state):
    
    import math
    index = index - 1
    
    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        S = 0.37464 + 1.54226 * w[i] - 0.26992 * (w[i]) ** (2)
        Tr = T / Tc[i]
        k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
        Ac[i] = 0.45724 * (R * Tc[i]) ** (2) * k / pc[i]
        Bc[i] = 0.0778 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)
    Ai = Ac[index] * pressure / (R * T) ** 2
    Bi = Bc[index] * pressure / (R * T)

    pig = 3.14159265353589
	
    beta = - 1 + b
    gamma = a - 2*b - 3*(b) ** (2)
    delta = -a * b + (b)**(2) + (b)**(3)
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus
    
    aus = Bi * (zeta - 1) / b + a * (Bi / b - 2 * (Ai / a) ** (0.5)) * math.log((zeta + b * (1 + (2) ** (0.5))) / (zeta + b * (1 - (2) ** (0.5)))) / ((2) ** (1.5) * b) - math.log(zeta - b) 
    
    phi_PRmix = math.exp(aus)

    return phi_PRmix
    
    
    
def hR_VdWmix(Tc,pc,T,pressure,x,state):
    
    import math

    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        #Tr = T / Tc[i]
        Ac[i] = 0.421875 * (R * Tc[i]) ** (2) / pc[i]
        Bc[i] = 0.125 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)
    
    pig = 3.14159265353589
	
    beta = - 1 - b
    gamma = a
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus
    
    hR_VdWmix = (zeta - 1 - a / zeta)*R*T

    return hR_VdWmix

def hR_RKmix(Tc,pc,T,pressure,x,state):
    
    import math

    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        Tr = T / Tc[i]
        k = 1 / (Tr)**(0.5)
        Ac[i] = 0.42748 * (R * Tc[i]) ** (2) * k / pc[i]
        Bc[i] = 0.08664 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)
    
    pig = 3.14159265353589
	
    beta = -1
    gamma = a - b - (b) ** (2)
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus
    
    hR_RKmix = (zeta - 1 - 1.5 * a / b * math.log((zeta + b)/zeta))*R*T

    return hR_RKmix
 
def hR_RKSmix(Tc,pc,w,T,pressure,x,state):
    
    import math

    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    S = [0 for i in range(len(Tc))]
    Tr = [0 for i in range(len(Tc))]
    k = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        S[i] = 0.48 + 1.574 * w[i] - 0.176 * (w[i]) ** (2)
        Tr[i] = T / Tc[i]
        k[i] = (1 + S[i] * (1 - (Tr[i]) ** (0.5))) ** (2)
        Ac[i] = 0.42748 * (R * Tc[i]) ** (2) * k[i] / pc[i]
        Bc[i] = 0.08664 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)
    
    pig = 3.14159265353589
	
    beta = -1
    gamma = a - b - (b) ** (2)
    delta = -a * b
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus
    
    e = 0
    
    for i in range(len(Tc)):
         e = e + x[i] * S[i] * (Ac[i] * Tr[i] / k[i])**(0.5)

    e = e / (amix)**(0.5)
    
    hR_RKSmix = (zeta - 1 - ( 1 + e ) * a / b * math.log((zeta + b)/zeta))*R*T

    return hR_RKSmix
 
def hR_PRmix(Tc,pc,w,T,pressure,x,state):
    
    import math

    amix = 0
    bmix = 0
    R = 8.3144621
    
    Ac = [0 for i in range(len(Tc))]
    Bc = [0 for i in range(len(Tc))]
    S = [0 for i in range(len(Tc))]
    Tr = [0 for i in range(len(Tc))]
    k = [0 for i in range(len(Tc))]
    
    for i in range(len(Tc)):
        S[i] = 0.37464 + 1.54226 * w[i] - 0.26992 * (w[i]) ** (2)
        Tr[i] = T / Tc[i]
        k[i] = (1 + S[i] * (1 - (Tr[i]) ** (0.5))) ** (2)
        Ac[i] = 0.45724 * (R * Tc[i]) ** (2) * k[i] / pc[i]
        Bc[i] = 0.0778 * R * Tc[i] / pc[i]

    for i in range(len(Tc)):
        for j in range(len(Tc)):
            amix = amix + x[i] * x[j] * (Ac[i] * Ac[j]) ** (0.5)
            bmix = bmix + x[i] * x[j] * (Bc[i] + Bc[j]) / 2

    a = amix * pressure / (R * T) ** (2)
    b = bmix * pressure / (R * T)
    
    pig = 3.14159265353589
	
    beta = - 1 + b
    gamma = a - 2*b - 3*(b) ** (2)
    delta = -a * b + (b)**(2) + (b)**(3)
     
    p = gamma - (beta) ** (2) / 3
    q = 2 * (beta) ** (3) / 27 - beta * gamma / 3 + delta
    det = (q) ** (2) / 4 + (p) ** (3) / 27
   
    if det >= 0:
        zeta1 = -q / 2 + (det) ** (0.5)
        zeta2 = -q / 2 - (det) ** (0.5)

        coeff1 = 1
        coeff2 = 1
        
        if zeta1 < 0:
            zeta1 = -zeta1
            coeff1 = -1

        if zeta2 < 0:
            zeta2 = -zeta2
            coeff2 = -1

        u = coeff1 * (zeta1) ** (1 / 3)
        v = coeff2 * (zeta2) ** (1 / 3)
        y1 = u + v
        aus = y1 - beta / 3
        
    else:
        if q < 0:
            teta = math.atan(-2 * (-det) ** (0.5) / q)
        else:
            teta = pig + math.atan(-2 * (-det) ** (0.5) / q)
      
        y1 = 2 * (-p / 3) ** (0.5) * math.cos(teta / 3)
        y2 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 2 * pig) / 3)
        y3 = 2 * (-p / 3) ** (0.5) * math.cos((teta + 4 * pig) / 3)
        zeta1 = y1 - beta / 3
        zeta2 = y2 - beta / 3
        zeta3 = y3 - beta / 3
		
        aus = zeta1
        
        if state == "L":
            #LIQUID PHASE
            if zeta2 < aus:
                aus = zeta2
            if zeta3 < aus:
                aus = zeta3
        elif state == "V":
            #VAPOUR PHASE
            if zeta2 > aus:
                aus = zeta2
            if zeta3 > aus:
                aus = zeta3

    zeta = aus
    
    e = 0
    
    for i in range(len(Tc)):
         e = e + x[i] * S[i] * (Ac[i] * Tr[i] / k[i])**(0.5)

    e = e / (amix)**(0.5)
    
    hR_PRmix = (zeta - 1 - ( 1 + e ) * a / b / (2)**(1.5) * math.log((zeta + b * (1 + (2)**(0.5)))/(zeta + b * (1 - (2)**(0.5)))))*R*T

    return hR_PRmix
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
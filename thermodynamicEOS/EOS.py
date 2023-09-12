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

class EOS:

    # checking arguments
    def check(state,p,T,pc,Tc,x=0,w=0):

        fatalError = 0
        p = EOS.checkPressure(p)
        T = EOS.checkTemperature(T)
        state = EOS.checkState(state)
        x = EOS.checkMoleFractions(x)
        index = EOS.checkIndex(len(Tc),index=1)

        if not(isinstance(x,(int,float))):
            if w==0:
                if not(len(x)==len(pc) and len(pc)==len(Tc)):
                    fatalError = 100
                    print("\n-----------------------------------------------")
                    print("\t\t!FATAL ERROR!")
                    print("  Tc,pc,x arrays must be of the same size")
                    print("-----------------------------------------------")
            else:
                if not(len(x)==len(pc) and len(pc)==len(Tc) and len(x)==len(w)):
                    fatalError = 101
                    print("\n-----------------------------------------------")
                    print("\t\t!FATAL ERROR!")
                    print("  Tc,pc,w,x arrays must be of the same size")
                    print("-----------------------------------------------")

        return p,T,x,state,index,fatalError
    
    def checkPressure(p):
        if p<=0:
            print("\n-----------------------------------------------")
            print("\t\t!ERROR!")
            print("  pressure should be greater than 0")
            print("  pressure set to 1 bar (1e5 Pa) by default")
            print("-----------------------------------------------")
            return 1e5
        else:
            return p
    def checkTemperature(T):
        if T<=0:
            print("\n-----------------------------------------------")
            print("\t\t!ERROR!")
            print("  temperature should be greater than 0")
            print("  temperature set to 300 K by default")
            print("-----------------------------------------------")
            return 300.
        else:
            return T
    def checkState(state):
        if state=="V" or state=="L":
            return state
        else:
            print("\n-----------------------------------------------")
            print("\t\t!ERROR!")
            print("  state should be either 'V' or 'L'")
            print("  state defined as 'V' by default")
            print("-----------------------------------------------")
            return "V"       
    def checkMoleFractions(x):
        x = abs(x)
        if sum(x)==1.:
            print("\n-----------------------------------------------")
            print("\t\t!ERROR!")
            print("  x elements should be positive and add up to 1")
            print("  x vector normalized by default")
            print("-----------------------------------------------")
            return x
        else:
            return x/sum(x)
    def checkIndex(N,index):
        index = int(index)
        if index<=0 or index>N:
            print("\n-----------------------------------------------")
            print("\t\t!ERROR!")
            print("  index should be an integer between 1 and {:d}".format(N))
            print("  pressure set to 1 bar (1e5 Pa) by default")
            print("-----------------------------------------------\n")
            index = 1
        return index
    


    def zetaCubicSolver(beta,gamma,delta,state):
        import math
        pig = 3.14159265353589
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
        return aus
    
    
                                                                              
    def zeta_VdW(Tc,pc,T,pressure,state):
    
        R = 8.3144621
        
        a_min = 0.421875 * (R * Tc) ** (2) / pc
        b_min = 0.125 * R * Tc / pc
        a = a_min * pressure / (R * T) ** (2)
        b = b_min * pressure / (R * T)   
        
        beta = - 1 - b
        gamma = a
        delta = -a * b
             
        zeta_VdW = EOS.zetaCubicSolver(beta,gamma,delta,state)
        return zeta_VdW
        
    def zeta_RK(Tc,pc,T,pressure,state):
     
        R = 8.3144621
        
        Tr = T / Tc
        k = 1/(Tr) ** (0.5)
        a_min = 0.42748 * (R * Tc) ** (2) * k / pc
        b_min = 0.08664 * R * Tc / pc
        a = a_min * pressure / (R * T) ** (2)
        b = b_min * pressure / (R * T)
        
        beta = - 1
        gamma =  a - b - (b) ** (2)
        delta = -a * b
         
        zeta_RK = EOS.zetaCubicSolver(beta,gamma,delta,state)
        return zeta_RK
    
    def zeta_RKS(Tc,pc,w,T,pressure,state):
    
        R = 8.3144621
        
        Tr = T / Tc
        S = 0.48 + 1.574 * w - 0.176 * (w) ** (2)
        k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
        a_min = 0.42748 * (R * Tc) ** (2) * k / pc
        b_min = 0.08664 * R * Tc / pc
        a = a_min * pressure / (R * T) ** (2)
        b = b_min * pressure / (R * T)
        
        beta = - 1
        gamma = a - b - (b) ** (2)
        delta = - a * b
         
        zeta_RKS = EOS.zetaCubicSolver(beta,gamma,delta,state)
        return zeta_RKS
        
    def zeta_PR(Tc,pc,w,T,pressure,state):
    
        R = 8.3144621
        
        Tr = T / Tc
        S = 0.37464 + 1.54226 * w - 0.26992 * w ** (2)
        k = (1 + S * (1 - (Tr) ** (0.5))) ** (2)
        a_min = 0.45724 * (R * Tc) ** (2) * k / pc
        b_min = 0.0778 * R * Tc / pc
        a = a_min * pressure / (R * T) ** (2)
        b = b_min * pressure / (R * T)
        
        beta = - 1 + b
        gamma = a - 2 * b - 3 * (b) ** (2)
        delta = -a * b + (b) ** (2) + (b) ** (3)
         
        zeta_PR = EOS.zetaCubicSolver(beta,gamma,delta,state)
        return zeta_PR
    
    
    
    def phi_VdW(Tc,pc,T,pressure,state):
        import math
        R = 8.3144621
        
        #Tr = T / Tc
        a_min = 0.421875 * (R * Tc) ** (2) / pc
        b_min = 0.125 * R * Tc / pc
        a = a_min * pressure / (R * T) ** (2)
        b = b_min * pressure / (R * T)
        
        beta = - 1 - b
        gamma = a
        delta = - a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
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
        
        beta = - 1
        gamma = a - b - (b) ** (2)
        delta = - a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
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
        
        beta = - 1
        gamma = a - b - (b) ** (2)
        delta = - a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
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
        
        beta = - 1 + b
        gamma = a - 2 * b - 3 * (b) ** (2)
        delta = -a * b + (b) ** (2) + (b) ** (3)
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        aus = zeta - 1 - a * math.log((zeta + b) / zeta) / b - math.log(zeta - b)
        phi_PR = math.exp(aus)
        return phi_PR
        
    
    
    def hR_VdW(Tc,pc,T,pressure,state):
        


        R = 8.3144621
        
        a_min = 0.421875 * (R * Tc) ** (2) / pc
        b_min = 0.125 * R * Tc / pc
        a = a_min * pressure / (R * T) ** (2)
        b = b_min * pressure / (R * T)
        
        beta = - 1 - b
        gamma = a
        delta = - a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
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
     
        beta = - 1
        gamma = a - b - (b) ** (2)
        delta = - a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
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
     
        beta = - 1
        gamma = a - b - (b) ** (2)
        delta = - a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
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
    
        beta = - 1 + b
        gamma = a - 2 * b - 3 * (b) ** (2)
        delta = - a * b + (b) ** (2) + (b) ** (3)
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
        e = S * (Tr / k)**(0.5)
        hR_PR = (zeta - 1 - ( 1 + e ) * a / b / (2)**(1.5) * math.log((zeta + b * (1 + (2)**(0.5)))/(zeta + b * (1 - (2)**(0.5)))))*R*T
        return hR_PR
        
      
      
    def zeta_VdWmix(Tc,pc,T,pressure,x,state):
        
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
      
        beta = - 1 - b
        gamma = a
        delta = -a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
    
        zeta_VdWmix = zeta
        return zeta_VdWmix
    
    def zeta_RKmix(Tc,pc,T,pressure,x,state):
        
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
     
        beta = -1
        gamma = a - b - (b) ** (2)
        delta = -a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
    
        zeta_RKmix = zeta
        return zeta_RKmix
    
    def zeta_RKSmix(Tc,pc,w,T,pressure,x,state):     
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
    
        beta = -1
        gamma = a - b - (b) ** (2)
        delta = -a * b
    
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
    
        zeta_RKSmix = zeta
        return zeta_RKSmix
    
    def zeta_PRmix(Tc,pc,w,T,pressure,x,state):
    
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
    
        beta = - 1 + b
        gamma = a - 2*b - 3*(b) ** (2)
        delta = -a * b + (b)**(2) + (b)**(3)
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
    
        zeta_PRmix = zeta
        return zeta_PRmix
    
    
    
    def phi_VdWmix(Tc,pc,index,T,pressure,x,state):
        import math

        state = EOS.checkState(state)

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
    
        beta = - 1 - b
        gamma = a
        delta = -a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
        aus = Bi/(zeta-b)-2*(Ai*a)**(0.5)/zeta-math.log(zeta-b) 
        
        phi_VdWmix = math.exp(aus)
    
        return phi_VdWmix
    
    def phi_RKmix(Tc,pc,index,T,pressure,x,state):
        import math

        state = EOS.checkState(state)

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
     
        beta = -1
        gamma = a - b - (b) ** (2)
        delta = -a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
        aus = Bi * (zeta - 1) / b + a * (Bi / b - 2 * (Ai / a) ** (0.5)) * math.log((zeta + b) / zeta) / b - math.log(zeta - b) 
        
        phi_RKmix = math.exp(aus)
    
        return phi_RKmix
        
    def phi_RKSmix(Tc,pc,w,index,T,pressure,x,state):
        import math

        state = EOS.checkState(state)

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
    
        beta = -1
        gamma = a - b - (b) ** (2)
        delta = -a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
        aus = Bi * (zeta - 1) / b + a * (Bi / b - 2 * (Ai / a) ** (0.5)) * math.log((zeta + b) / zeta) / b - math.log(zeta - b) 
        
        phi_RKSmix = math.exp(aus)
    
        return phi_RKSmix
        
    def phi_PRmix(Tc,pc,w,index,T,pressure,x,state):
        import math

        state = EOS.checkState(state)

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
    
        beta = - 1 + b
        gamma = a - 2*b - 3*(b) ** (2)
        delta = -a * b + (b)**(2) + (b)**(3)
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
        aus = Bi * (zeta - 1) / b + a * (Bi / b - 2 * (Ai / a) ** (0.5)) * math.log((zeta + b * (1 + (2) ** (0.5))) / (zeta + b * (1 - (2) ** (0.5)))) / ((2) ** (1.5) * b) - math.log(zeta - b) 
        
        phi_PRmix = math.exp(aus)
    
        return phi_PRmix
        
        
        
    def hR_VdWmix(Tc,pc,T,pressure,x,state):
        
        state = EOS.checkState(state)

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
        
        beta = - 1 - b
        gamma = a
        delta = -a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
        hR_VdWmix = (zeta - 1 - a / zeta)*R*T
    
        return hR_VdWmix
    
    def hR_RKmix(Tc,pc,T,pressure,x,state):
        import math

        state = EOS.checkState(state)
    
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
    
        beta = -1
        gamma = a - b - (b) ** (2)
        delta = -a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
        hR_RKmix = (zeta - 1 - 1.5 * a / b * math.log((zeta + b)/zeta))*R*T
    
        return hR_RKmix
     
    def hR_RKSmix(Tc,pc,w,T,pressure,x,state):
        import math

        state = EOS.checkState(state)
    
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
        
        beta = -1
        gamma = a - b - (b) ** (2)
        delta = -a * b
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
        e = 0
        
        for i in range(len(Tc)):
             e = e + x[i] * S[i] * (Ac[i] * Tr[i] / k[i])**(0.5)
    
        e = e / (amix)**(0.5)
        
        hR_RKSmix = (zeta - 1 - ( 1 + e ) * a / b * math.log((zeta + b)/zeta))*R*T
    
        return hR_RKSmix
     
    def hR_PRmix(Tc,pc,w,T,pressure,x,state):
        import math

        state = EOS.checkState(state)
    
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
        
        beta = - 1 + b
        gamma = a - 2*b - 3*(b) ** (2)
        delta = -a * b + (b)**(2) + (b)**(3)
         
        zeta = EOS.zetaCubicSolver(beta,gamma,delta,state)
        
        e = 0
        
        for i in range(len(Tc)):
             e = e + x[i] * S[i] * (Ac[i] * Tr[i] / k[i])**(0.5)
    
        e = e / (amix)**(0.5)
        
        hR_PRmix = (zeta - 1 - ( 1 + e ) * a / b / (2)**(1.5) * math.log((zeta + b * (1 + (2)**(0.5)))/(zeta + b * (1 - (2)**(0.5)))))*R*T
    
        return hR_PRmix        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
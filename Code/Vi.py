import math

def GetV_i(index, r0, r1, r2, u1, u2, N):
    rho3 = 1.0 - (r0*(1-u1) / (r2 + r1 * u2))
    rho3 /= 1.0 - math.pow(r0*(1.0-u1) / (r2 + r1 * u2), N)
    
    V = []
    V_new = []
    V_lin = []
    
    for i in range(0,N+1):
        V.append(1.0 - float(i) / (N))
        V_new.append(V[i])
        V_lin.append(V[i])
        #print("V{0} = {1}".format(i, V[i]))
    
    MAX_ITER = 1000
    MIN_CHANGE = 1e-10
    
    V1_old = V[1]
    
    for it in range(0, MAX_ITER):
        for i in range(1, N):
            nxt = V[i+1]
            prv = V[i-1]
            P_i = float(i*(N-i)) / (i*r1*(1.0-u2) + r0*(N-i))
            new = (r1*(1.0-u2) * nxt + r0 * prv)
            new = new / ((r1*u2*rho3)/P_i + r1*(1.0-u2) + r0)
            V_new[i] = new
        
        for i in range(1,N):
            V[i] = V_new[i]        
        
        if abs(V[1] - V1_old) < MIN_CHANGE:
            #print("Converged in {0} steps.".format(i))
            break
        V1_old = V[1]
        
    return V[index]


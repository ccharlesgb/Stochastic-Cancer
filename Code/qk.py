import scipy.integrate

def GetLambdaK(k, N, r1, r2, u2):
    res = (N - k) * ((N-k)*r1*u2 + r2*k) / ((N-k)*r1 + r2*k)
    
    return res


def GetMuK(k, N, r1, r2, u2):
    res = k * ((N-k)*r1*(1-u2)) / ((N-k)*r1 + r2*k)
    
    return res


def GetQk(t, r0, r1, r2, u1, u2, N):
    q = []
    q_lin = []
    
    for i in range(0,N+1):
        q.append(float(i) / N)
        q_lin.append(q[i])
    
    
    MAX_ITER = 1000
    MIN_CHANGE = 1e-5
    
    V1_old = V[1]
    
    for it in range(0, MAX_ITER):
        for i in range(1, N):
            nextq = q[i + 1]
            prevq = q[i - 1]
            curq = q[i]
            curlambda = GetLambdaK(k, N, r1, r2, u2)
            curmu = GetMuK(k, N, r1, r2, u2)
            
            qkdot = curlambda * nextq - (curlambda + curmu)*curq + curmu * prevq
            newQ = qkdot * t
            q[i] = newQ
    
    return q[0]
            

'''
for it in range(0, MAX_ITER):
    for i in range(1, N-1):
        nxt = V[i+1]
        prv = V[i-1]
        P_i = float(i*(N-i)) / (i*r1*(1.0-u2) + r0*(N-i))
        new = (r1*(1-u2) * nxt + r0 * prv)
        new = new / ((r1*u2*rho3)/P_i + r1*(1-u2) + r0)
        V[i] = new
    
    if abs(V[1] - V1_old) < MIN_CHANGE:
        #print("Converged in {0} steps.".format(i))
        break
    V1_old = V[1]
    
return V[index]
'''

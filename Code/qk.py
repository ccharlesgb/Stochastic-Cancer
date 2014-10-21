import scipy.integrate
import matplotlib.pyplot as plt
import math

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
        q.append((float(i) / N))
        q_lin.append(q[i])
    
    
    delta_t = 0.1
    MAX_ITER = int(t / delta_t)
    MIN_CHANGE = 1e-5
    
    curt = 0.0
    for it in range(0, MAX_ITER):
        for k in range(0, N):
            nextq = q[k + 1]
            if k > 0:
                prevq = q[k - 1]
            else:
                prevq = 0.0
            curq = q[k] 
            curlambda = GetLambdaK(k, N, r1, r2, u2)
            if k > 0:
                curmu = GetMuK(k, N, r1, r2, u2)
            else:
                curmu = 0.0
            
            qkdot = curlambda * nextq - (curlambda + curmu)*curq + curmu * prevq
            newQ = q[k] + qkdot * delta_t
            q[k] = newQ
        
        curt = it * delta_t
    
    return q[0]
           
datapointcount = 25           
           
dataPointsq=[]
dataPointst=[]
dataPointsG = []

t = 100.0
r0 = 1.0
r1 = 1.01
r2 = 1.0
u1 = 0.1
u2 = 0.1

N = 10

for i in range(0, datapointcount):
    t = float(i) / (datapointcount - 1) * 100.0
    q = GetQk(t, r0, r1, r2 ,u1 ,u2, N)
    dataPointsq.append(q)
    dataPointst.append(t)
    
    c=N*u2*(1.0-r1/r2)/(1.0-(math.pow((r1/r2) , N)))
    curG = 1.0 - math.exp(-c*t)
    dataPointsG.append(curG)
    
plt.plot(dataPointst, dataPointsq)
plt.plot(dataPointst, dataPointsG)
plt.show()
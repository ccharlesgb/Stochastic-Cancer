import scipy.integrate
import matplotlib.pyplot as plt
import math

def GetLambdaK(k, N, r1, r2, u2):
    res = (N - k) * ((N-k)*r1*u2 + r2*k) / ((N-k)*r1 + r2*k)
    return res


def GetMuK(k, N, r1, r2, u2):
    res = k * ((N-k)*r1*(1-u2)) / ((N-k)*r1 + r2*k)
    return res


def GetQ0(t, r0, r1, r2, u1, u2, N):
    return GetQk(0, t, r0, r1, r2, u1, u2, N)

def GetQk(k_return, t, r0, r1, r2, u1, u2, N):
    q = []
    q_dots = []
    q_lin = []
    
    for i in range(0,N+1):
        q.append(0.0)
        q_dots.append(0.0)
        q_lin.append(q[i])
    
    q[N] = 1.0
    delta_t = 0.1
    MAX_ITER = int(t / delta_t)

    curt = 0.0
    for it in range(0, MAX_ITER + 1):
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
            
            q_dots[k] = curlambda * nextq - (curlambda + curmu)*curq + curmu * prevq
        
        if (it == MAX_ITER):
            delta_t = t - MAX_ITER * delta_t
            #print("Final Step: {0}".format(delta_t))
        curt = curt + delta_t
        for k in range(0,N):
            newQ = q[k] + q_dots[k] * delta_t
            q[k] = newQ
    
   # print("Final T {0}   {1}".format(curt, t))
    return q[k_return]
    
'''
datapointcount = 100           
           
dataPointsq=[]
dataPointst=[]
dataPointsG = []

dataPointsk = []
dataPointsQk = []

t = 100.0
r0 = 1.0
r1 = 0.5
r2 = 1.0
u1 = 0.1
u2 = 0.1

N = 10


for k in range(0, N + 1):
    dataPointsk.append(k)
    dataPointsQk.append(GetQk(k, t, r0, r1, r2, u1, u2, N))
    
plt.figure()
plt.plot(dataPointsk, dataPointsQk)

plt.figure()

for i in range(0, datapointcount):
    t = float(i) / (datapointcount - 1) * 100.0
    q = GetQk(0, t, r0, r1, r2 ,u1 ,u2, N)
    dataPointsq.append(q)
    dataPointst.append(t)
    
    c=N*u2*(1.0-r1/r2)/(1.0-(math.pow((r1/r2) , N)))
    curG = 1.0 - math.exp(-c*t)
    dataPointsG.append(curG)
    
plt.plot(dataPointst, dataPointsq, label = "qk")
plt.plot(dataPointst, dataPointsG, label = "G(t)", ls = '--')
plt.legend()
plt.show()
'''

r0 = 1.0
r1 = 1.0
r2 = 1.0

u1 = 0.1
u2 = 0.01

N = 100

rho3 = 1.0 - (r0*(1-u1) / (r2 + r1 * u2))
rho3 /= 1.0 - pow(r0*(1-u1) / (r2 + r1 * u2), N)

V = []

for i in range(0,N):
    V.append(1.0 - float(i) / N)

MAX_ITER = 100

for it in range(0, MAX_ITER):
    for i in range(1, N-1):
        nxt = V[i+1]
        prv = V[i-1]
        P_i = float(i*(N-i)) / (i*r1*(1.0-u2) + r0*(N-i))
        new = (r1*(1-u2) * nxt + r0 * prv)
        new = new / ((r1*u2*rho3)/P_i + r1*(1-u2) + r0)
        V[i] = new
    print("New V1 = {0}".format(V[1]))
    

   
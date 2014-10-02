import matplotlib.pyplot as plt
import SimTools

mySim = SimTools.Gillespie(100)
mySim.timeLimit = 100
mySim.u2 = 0.01
dataPointCount = 25

simsPerDataPoint = 1000

dataPointsX = []
dataPointsY = []
for i in range(0, dataPointCount):
    dataPointsX.insert(i,0.0)
    dataPointsY.insert(i,0.0)

for curPoint in range(0, dataPointCount):
    fixationCounts = 0
    mySim.r1 = 0.2 + (3.0-0.2) * float(curPoint)/float(dataPointCount) #sweep type-1 fitness from 0.2 to 3.0
    print("Type-1 Fitness: %f" % mySim.r1)
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()

        if mySim.n2 == mySim.N:
            fixationCounts += 1
    dataPointsX[curPoint] = mySim.r1
    dataPointsY[curPoint] = float(fixationCounts) / float(simsPerDataPoint)

for i in range(0, dataPointCount):
    print("r1: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))


plt.plot(dataPointsX, dataPointsY)
plt.xlabel("Type 1 Fitness r1")
plt.ylabel("Type 2 Fixation %")
plt.show()

   
filename = "sim_N={0}_r0={1}_r1={2}_r2={3}_u1={4}_u2={5}_SPDP={6}".format(mySim.N, mySim.r0, mySim.r1, mySim.r2, mySim.u1, mySim.u2, simsPerDataPoint)

f = open(filename, 'w')

for i in range(0, dataPointCount):
    f.write(str(dataPointsX[i]) + "    " + str(dataPointsY[i]) + "\n")

f.close()



        

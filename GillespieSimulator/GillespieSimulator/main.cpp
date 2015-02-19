#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "Gillespie.h"

int main()
{
	Gillespie myGill = Gillespie();
	int TYPE_COUNT = 20;
	GillespieParam myParam = GillespieParam(TYPE_COUNT);
	myParam.n0[0] = 10;
	myGill.Hook(&myParam);

	GillespieHist myHist(TYPE_COUNT);
	myGill.SetHistory(&myHist);

	myGill.SetTimeLimit(10000.0f);

	//myGillespie.RECORD_TAU_INFO = 1

	myParam.n0[0] = 1e6f;

	for (int i = 0; i < TYPE_COUNT; i++)
	{
		myParam.r[i] = std::pow(1.0 + 0.01, i);
		myParam.u[i] = 1.0 / myParam.n0[0];
	}
    
	myGill.SetEpsilon(0.3f);

	myGill.Simulate();


	for (int i=0; i<myGill.GetSimSteps(); i+= 10)
	{
		std::cout << myHist.tHist[i] << "  ";
		for (int pop = 0; pop < TYPE_COUNT; pop++)
		{
			std::cout << myHist.nHist[pop][i] << "  ";
		}
		std::cout << "\n";
	}
		
	float test = 0.0f;
	std::cin >> test;
	return 0;
}
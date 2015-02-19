#include "GillespieParam.h"
#include "Gillespie.h"

void GillespieParam::Hook(Gillespie *gill)
{
	for (int i=0; i<mTypeCount; i++)
	{
		for (int j=0; j<mTypeCount; j++)
		{
			if (i == j)
				continue;
			std::cout << "Adding event " << i << " " << j << "\n";
			gill->AddCallback(EventTIJ(i,j));
		}
	}
}

std::vector<int>* GillespieParam::EventTIJ(int i, int j)
{
	std::vector<int> *arr = new std::vector<int>();
    for (int x = 0; x < mTypeCount; x++)
	{
        arr->push_back(0);
        if (x == i)
            (*arr)[x] = -1;
        else if (x == j)
            (*arr)[x] = 1;
	}
    return arr;
}

float GillespieParam::GetTIJ(int i, int j)
{
	if (n[i] == 0) //Quick get out case
		return 0.0;
	float top = 0.0;
	if (j == 0)
		top = n[i] * (r[j] * (1 - u[j]) * n[j]);
	else if (j == mTypeCount - 1)
		top = n[i] * (r[j]*n[j] + r[j-1] * u[j-1] * n[j-1]);
	else
		top = n[i] * (r[j]*(1.0 - u[j])*n[j] + r[j-1]*u[j-1]*n[j-1]);
        
	return top / mAvgFit;
}

float GillespieParam::GetAvgFit()
{
	float tot = 0.0;
    for (int i = 0; i < mTypeCount; i++)
        tot += r[i] * n[i];
    return tot;
}
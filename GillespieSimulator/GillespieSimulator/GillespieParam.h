
#pragma once

#include <iostream>
#include <vector>

class Gillespie;

class GillespieParam
{
public:
	int mTypeCount;
	std::vector<int> n;
	std::vector<int> n0;
	std::vector<float> r;
	std::vector<float> u;
	int mN;
	float mAvgFit;

	GillespieParam(int typeCount)
	{
		mTypeCount = typeCount;
		for (int i = 0; i < typeCount; i++)
		{
			n.push_back(0);
			n0.push_back(0);
			r.push_back(1.0f);
			u.push_back(0.01f);
		}
	};

	~GillespieParam() {};

	void Reset()
	{
		mN = 0;
		for (int i = 0; i < mTypeCount; i++)
		{
			n[i] = n0[i];
			mN += n0[i];
		}
	}

	virtual void Hook(Gillespie *gill);

	virtual void PreSim() {mAvgFit = GetAvgFit();};
	virtual void PostSim() {};

	std::vector<int>* EventTIJ(int i, int j);
	float GetTIJ(int i, int j);
	float GetAvgFit();
	
};

class GillespieHist
{
public:
	int mTypeCount;
	std::vector<float> tHist;
	std::vector<std::vector<int>> nHist;


    GillespieHist(int typeCount)
	{
        mTypeCount = typeCount;
		for (int i=0; i<typeCount; i++)
		{
			nHist.push_back(*(new std::vector<int>));
		}
        ClearFrames();
	}

	~GillespieHist()
	{

	}
      
	void Reserve(int size)
	{
		tHist.reserve(size);
		for (int i=0; i<mTypeCount; i++)
		{
			nHist[i].reserve(size);
		}
	}

    void ClearFrames()
	{
        tHist.clear();
		for (int i=0; i<nHist.size(); i++)
		{
			nHist[i].clear();
		}
	}

    void RecordFrame(float time, GillespieParam *param)
	{
        tHist.push_back(time);
        for (int i=0; i<mTypeCount; i++)
		{
            nHist[i].push_back(param->n[i]);
		}
	}
         
    //def GetDictionary(self):
        //runDict = dict()
        //return runDict
};
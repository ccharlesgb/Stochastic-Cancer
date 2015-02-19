
#pragma once

#include <iostream>
#include <vector>
#include <random>
#include "GillespieParam.h"

class GillespieHist;

typedef float (*RateCallback)();

//Gillespie Simulator
class Gillespie
{
private:
	std::vector<RateCallback> rateCallbacks;
	std::vector<std::vector<int>> stateChanges;
	std::vector<float> rateCache;
	int rateCallbackCount;

	std::default_random_engine mRandGen;
	std::poisson_distribution<int> mPoission;

	std::vector<int> eventCount;
	
	GillespieParam *mParams;
	GillespieHist *mHist;

	float mLambda;
	int mSimSteps;
	float mCurTime;
	float mTimeLimit;
	float mEpsilon;

	float mTau;
	std::vector<float> tauCache;

	virtual void PreSim(GillespieParam *param) {};
	virtual void PostSim(GillespieParam *param) {};


public:
	Gillespie();
	~Gillespie();

	void Reset();

	void AddCallback(std::vector<int> *stateChange);

	float GetAuxMu(int popID);
	float GetAuxSigma(int popID);

	void GetTimeStep();

	void Hook(GillespieParam* param)
	{
		mParams = param;
		for (int i = 0; i < param->n.size(); i++)
		{
			tauCache.push_back(0.0);
		}
		param->Hook(this);
	};

	void UpdateRates();

	void DoEvents();

	void Simulate();

	void SetHistory(GillespieHist *hist)
	{
		mHist = hist;
	}

	int GetSimSteps() {return mSimSteps;};

	void SetEpsilon(float ep) {mEpsilon = ep;};
	void SetTimeLimit(float tl) {mTimeLimit = tl;};

};
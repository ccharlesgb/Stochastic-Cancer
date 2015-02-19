#include "Gillespie.h"
#include <math.h>
#include <algorithm>

Gillespie::Gillespie()
{
	rateCallbackCount = 0;
	mParams = NULL;
	mLambda = 0.0;
	mSimSteps = 0;
	mCurTime = 0;
	mTimeLimit = 100.0;
	mEpsilon = 0.05;
	mHist = NULL;
	mTau = -1.0;
}

Gillespie::~Gillespie()
{

}

void Gillespie::Reset()
{
	mSimSteps = 0;
	mCurTime = 0.0;
}

void Gillespie::AddCallback(std::vector<int> *stateChange)
{
	stateChanges.push_back(*stateChange);
	rateCallbackCount++;
	rateCache.push_back(0.0);
	eventCount.push_back(0.0);
}

float Gillespie::GetAuxMu(int popID)
{
	float sum = 0.0;
	for (int j = 0; j < rateCallbackCount; j++)
	{
		sum += stateChanges[j][popID] * rateCache[popID];
	}
	if (sum == 0)
	{
		sum = 1e-10f;
	}
	return sum;
}

float Gillespie::GetAuxSigma(int popID)
{
	float sum = 0.0;
	for (int j = 0; j < rateCallbackCount; j++)
	{
		sum += stateChanges[j][popID] * stateChanges[j][popID] * rateCache[j];
	}
	if (sum == 0)
		sum = 1e-10f;
	return sum;
}

void Gillespie::GetTimeStep()
{
	for (int i=0; i < mParams->mTypeCount; i++)
	{
		float comp = std::max<float>(mEpsilon*mParams->n[i], 1.0f);
		float tauElem = std::min<float>(comp / std::abs(GetAuxMu(i)) , (comp * comp) / GetAuxSigma(i) );
		tauCache[i] = tauElem;
	}
	mTau = *std::min_element(tauCache.begin(), tauCache.end());
	//mTau = 1.0f;
}

void Gillespie::UpdateRates()
{
    mLambda = 0.0;
    int index = 0;
	for (int i = 0; i < mParams->mTypeCount; i++)
	{
		for (int j = 0; j < mParams->mTypeCount; j++)
		{
            if (i == j)
                continue;
            rateCache[index] = mParams->GetTIJ(i,j);
            mLambda += rateCache[index];
                
            index += 1;
		}
	}

}

void Gillespie::DoEvents()
{
	bool goodFrame  = false;
	while (goodFrame == false)
	{
		goodFrame = true;
		for (int i = 0; i < rateCallbackCount; i++)
		{
			float newMean = rateCache[i] * mTau;
			if (newMean == 0)
				eventCount[i] = 0;
			else
			{
				mPoission.param(rateCache[i] * mTau);
				eventCount[i] = mPoission(mRandGen);
			}
			for (int pop = 0; pop < mParams->mTypeCount; pop++)
			{
				mParams->n[pop] += stateChanges[i][pop] * eventCount[i];
				if (mParams->n[pop] < 0)
				{
					goodFrame = false;
				}
			}
		}
		if (goodFrame == false)
		{
			//self.BAD_FRAME_COUNT += 1;
			for (int i = 0; i < rateCallbackCount; i++)
			{
				for (int pop = 0; pop < mParams->mTypeCount; pop++)
				{
					mParams->n[pop] -= stateChanges[i][pop] * eventCount[i];
				}
			}
			mTau /= 2.0;
		}
	}
}

void Gillespie::Simulate()
{
	mParams->Reset();
	Reset();

	if (mHist != NULL)
	{
		mHist->Reserve(mTimeLimit / mEpsilon);
		mHist->RecordFrame(mCurTime, mParams);
	}

	while (mCurTime < mTimeLimit)
	{
		mParams->PreSim();
		//TODO PRE SIM
		UpdateRates();
		if (mLambda == 0.0f)
		{
			mSimSteps++;
			if (mHist != NULL)
			{
				mHist->RecordFrame(mCurTime, mParams);
			}
			return;
		}
		GetTimeStep();
		DoEvents();
		mCurTime += mTau;
		mSimSteps++;
		//std::cout << mCurTime << "\n";
		mParams->PostSim();

		if (mSimSteps % 100 == 0)
			std::cout << mCurTime << "\n";

		if (mHist != NULL)
		{
			mHist->RecordFrame(mCurTime, mParams);
		}
	}
}
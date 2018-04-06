/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "UniformContactInhibitionCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "ApoptoticCellProperty.hpp"

UniformContactInhibitionCellCycleModel::UniformContactInhibitionCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mMinCellCycleDuration(12.0), // Hours
      mMaxCellCycleDuration(14.0),  // Hours
	  mQuiescentVolumeFraction(DOUBLE_UNSET),
      mEquilibriumVolume(DOUBLE_UNSET),
      mCurrentQuiescentOnsetTime(SimulationTime::Instance()->GetTime()),
      mCurrentQuiescentDuration(0.0)
{
}

UniformContactInhibitionCellCycleModel::UniformContactInhibitionCellCycleModel(const UniformContactInhibitionCellCycleModel& rModel)
   : AbstractSimpleCellCycleModel(rModel),
     mMinCellCycleDuration(rModel.mMinCellCycleDuration),
     mMaxCellCycleDuration(rModel.mMaxCellCycleDuration),
     mQuiescentVolumeFraction(rModel.mQuiescentVolumeFraction),
     mEquilibriumVolume(rModel.mEquilibriumVolume),
     mCurrentQuiescentOnsetTime(rModel.mCurrentQuiescentOnsetTime),
     mCurrentQuiescentDuration(rModel.mCurrentQuiescentDuration)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variable mCellCycleDuration is initialized in the
     * AbstractSimpleCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mCellCycleDuration is (re)set as soon as
     * InitialiseDaughterCell() is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* UniformContactInhibitionCellCycleModel::CreateCellCycleModel()
{
    return new UniformContactInhibitionCellCycleModel(*this);
}

void UniformContactInhibitionCellCycleModel::SetCellCycleDuration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCellCycleDuration = DBL_MAX;
    }
    else
    {
        mCellCycleDuration = mMinCellCycleDuration + (mMaxCellCycleDuration - mMinCellCycleDuration) * p_gen->ranf(); // U[MinCCD,MaxCCD]
    }
}

double UniformContactInhibitionCellCycleModel::GetMinCellCycleDuration()
{
    return mMinCellCycleDuration;
}

void UniformContactInhibitionCellCycleModel::SetMinCellCycleDuration(double minCellCycleDuration)
{
    mMinCellCycleDuration = minCellCycleDuration;
}

double UniformContactInhibitionCellCycleModel::GetMaxCellCycleDuration()
{
    return mMaxCellCycleDuration;
}

void UniformContactInhibitionCellCycleModel::SetMaxCellCycleDuration(double maxCellCycleDuration)
{
    mMaxCellCycleDuration = maxCellCycleDuration;
}

double UniformContactInhibitionCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

double UniformContactInhibitionCellCycleModel::GetAverageStemCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

void UniformContactInhibitionCellCycleModel::SetQuiescentVolumeFraction(double quiescentVolumeFraction)
{
    mQuiescentVolumeFraction = quiescentVolumeFraction;
}

double UniformContactInhibitionCellCycleModel::GetQuiescentVolumeFraction() const
{
    return mQuiescentVolumeFraction;
}

void UniformContactInhibitionCellCycleModel::SetEquilibriumVolume(double equilibriumVolume)
{
    mEquilibriumVolume = equilibriumVolume;
}

double UniformContactInhibitionCellCycleModel::GetEquilibriumVolume() const
{
    return mEquilibriumVolume;
}

double UniformContactInhibitionCellCycleModel::GetCurrentQuiescentDuration() const
{
    return mCurrentQuiescentDuration;
}

double UniformContactInhibitionCellCycleModel::GetCurrentQuiescentOnsetTime() const
{
    return mCurrentQuiescentOnsetTime;
}

void UniformContactInhibitionCellCycleModel::UpdateQuiescenceBasedOnCellVolume()
{
    assert(!(mpCell->HasCellProperty<ApoptoticCellProperty>()));
    assert(!mpCell->HasApoptosisBegun());

    if ((mQuiescentVolumeFraction == DOUBLE_UNSET) || (mEquilibriumVolume == DOUBLE_UNSET))
    {
        EXCEPTION("The member variables mQuiescentVolumeFraction and mEquilibriumVolume have not yet been set.");
    }

    // Get cell volume
    double cell_volume = mpCell->GetCellData()->GetItem("volume");

    // Freeze cell based on cell volume
    double dt = SimulationTime::Instance()->GetTimeStep();
    double quiescent_volume = mEquilibriumVolume * mQuiescentVolumeFraction;

    if (cell_volume < quiescent_volume)
    {
        // Update the duration of the current period of contact inhibition.
        mCurrentQuiescentDuration = SimulationTime::Instance()->GetTime() - mCurrentQuiescentOnsetTime;

        // Increase cell cycle duration by length of timestep, essentially freezing the cell cycle for this timestep
    	mCellCycleDuration = mCellCycleDuration + dt;

    }
    else
    {
        // Reset the cell's quiescent duration and update the time at which the onset of quiescent occurs
        mCurrentQuiescentDuration = 0.0;
        mCurrentQuiescentOnsetTime = SimulationTime::Instance()->GetTime();
    }

}

// Overridden ReadyToDivide
bool UniformContactInhibitionCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);
    UpdateQuiescenceBasedOnCellVolume();

    if (!mReadyToDivide)
    {
        if (GetAge() >= mCellCycleDuration )
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
}


void UniformContactInhibitionCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<QuiescentVolumeFraction>" << mQuiescentVolumeFraction << "</QuiescentVolumeFraction>\n";
    *rParamsFile << "\t\t\t<EquilibriumVolume>" << mEquilibriumVolume << "</EquilibriumVolume>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(UniformContactInhibitionCellCycleModel)

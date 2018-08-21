///*
//
//Copyright (c) 2005-2017, University of Oxford.
//All rights reserved.
//
//University of Oxford means the Chancellor, Masters and Scholars of the
//University of Oxford, having an administrative office at Wellington
//Square, Oxford OX1 2JD, UK.
//
//This file is part of Chaste.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of the University of Oxford nor the names of its
//   contributors may be used to endorse or promote products derived from this
//   software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//*/
//
//#include "UniformCellCycleModelWithQuiescence.hpp"
//#include "DifferentiatedCellProliferativeType.hpp"
//
//#include "ApoptoticCellProperty.hpp"
//
//UniformCellCycleModelWithQuiescence::UniformCellCycleModelWithQuiescence()
//    : AbstractSimpleCellCycleModel(),
//      mMinCellCycleDuration(12.0), // Hours
//      mMaxCellCycleDuration(14.0),  // Hours
//	  mCurrentHypoxicDuration(0.0),
//	  mHypoxicConcentration(0.4),
//	  mQuiescentConcentration(1.0),
//	  mCriticalHypoxicDuration(2.0)
//{
//	mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetTime();
//}
//
//UniformCellCycleModelWithQuiescence::UniformCellCycleModelWithQuiescence(const UniformCellCycleModelWithQuiescence& rModel)
//   : AbstractSimpleCellCycleModel(rModel),
//     mMinCellCycleDuration(rModel.mMinCellCycleDuration),
//     mMaxCellCycleDuration(rModel.mMaxCellCycleDuration),
//	 mCurrentHypoxicDuration(rModel.mCurrentHypoxicDuration),
//	 mCurrentHypoxiaOnsetTime(rModel.mCurrentHypoxiaOnsetTime),
//	 mHypoxicConcentration(rModel.mHypoxicConcentration),
//	 mQuiescentConcentration(rModel.mQuiescentConcentration),
//	 mCriticalHypoxicDuration(rModel.mCriticalHypoxicDuration)
//{
//    /*
//     * Initialize only those member variables defined in this class.
//     *
//     * The member variable mCellCycleDuration is initialized in the
//     * AbstractSimpleCellCycleModel constructor.
//     *
//     * The member variables mBirthTime, mReadyToDivide and mDimension
//     * are initialized in the AbstractCellCycleModel constructor.
//     *
//     * Note that mCellCycleDuration is (re)set as soon as
//     * InitialiseDaughterCell() is called on the new cell-cycle model.
//     */
//}
//
//AbstractCellCycleModel* UniformCellCycleModelWithQuiescence::CreateCellCycleModel()
//{
//    return new UniformCellCycleModelWithQuiescence(*this);
//}
//
//void UniformCellCycleModelWithQuiescence::SetCellCycleDuration()
//{
//    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
//
//    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
//    {
//        mCellCycleDuration = DBL_MAX;
//    }
//    else
//    {
//        mCellCycleDuration = mMinCellCycleDuration + (mMaxCellCycleDuration - mMinCellCycleDuration) * p_gen->ranf(); // U[MinCCD,MaxCCD]
//    }
//}
//
//double UniformCellCycleModelWithQuiescence::GetMinCellCycleDuration()
//{
//    return mMinCellCycleDuration;
//}
//
//void UniformCellCycleModelWithQuiescence::SetMinCellCycleDuration(double minCellCycleDuration)
//{
//    mMinCellCycleDuration = minCellCycleDuration;
//}
//
//double UniformCellCycleModelWithQuiescence::GetMaxCellCycleDuration()
//{
//    return mMaxCellCycleDuration;
//}
//
//void UniformCellCycleModelWithQuiescence::SetMaxCellCycleDuration(double maxCellCycleDuration)
//{
//    mMaxCellCycleDuration = maxCellCycleDuration;
//}
//
//double UniformCellCycleModelWithQuiescence::GetAverageTransitCellCycleTime()
//{
//    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
//}
//
//double UniformCellCycleModelWithQuiescence::GetAverageStemCellCycleTime()
//{
//    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
//}
//
//double UniformCellCycleModelWithQuiescence::GetCurrentHypoxicDuration() const
//{
//    return mCurrentHypoxicDuration;
//}
//
//void UniformCellCycleModelWithQuiescence::UpdateHypoxicDuration()
//{
//    assert(!(mpCell->HasCellProperty<ApoptoticCellProperty>()));
//    assert(!mpCell->HasApoptosisBegun());
//
//    // Get cell's oxygen concentration
//    double oxygen_concentration = mpCell->GetCellData()->GetItem("oxygen");
//
//    if (oxygen_concentration < mHypoxicConcentration)
//    {
//        // Update the duration of the current period of hypoxia
//        mCurrentHypoxicDuration = (SimulationTime::Instance()->GetTime() - mCurrentHypoxiaOnsetTime);
//
//        if (mCurrentHypoxicDuration > mCriticalHypoxicDuration)
//        {
//            /*
//             * This method is usually called within a CellBasedSimulation, after the CellPopulation
//             * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
//             * CellPropertyRegistry::Instance() here when adding the ApoptoticCellProperty, we would
//             * be creating a new CellPropertyRegistry. In this case the ApoptoticCellProperty cell
//             * count would be incorrect. We must therefore access the ApoptoticCellProperty via the
//             * cell's CellPropertyCollection.
//             */
//            boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
//                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
//            mpCell->AddCellProperty(p_apoptotic_property);
//        }
//    }
//    else
//    {
//        // Reset the cell's hypoxic duration and update the time at which the onset of hypoxia occurs
//        mCurrentHypoxicDuration = 0.0;
//        mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetTime();
//    }
//
//    // If cell is quiescent, freeze its cell cycle. We can do this by either decreasing its age or increasing its cell cycle duration. Here, we increase the cell cycle duration
//    if (oxygen_concentration < mQuiescentConcentration)
//    {
//    	// Increase cell cycle duration by length of timestep
//    	mCellCycleDuration = mCellCycleDuration + SimulationTime::Instance()->GetTimeStep();
//
//    }
//}
//
//double UniformCellCycleModelWithQuiescence::GetHypoxicConcentration() const
//{
//    return mHypoxicConcentration;
//}
//
//void UniformCellCycleModelWithQuiescence::SetHypoxicConcentration(double hypoxicConcentration)
//{
//    assert(hypoxicConcentration<=1.0);
//    assert(hypoxicConcentration>=0.0);
//    mHypoxicConcentration = hypoxicConcentration;
//}
//
//double UniformCellCycleModelWithQuiescence::GetQuiescentConcentration() const
//{
//    return mQuiescentConcentration;
//}
//
//void UniformCellCycleModelWithQuiescence::SetQuiescentConcentration(double quiescentConcentration)
//{
//    assert(quiescentConcentration <= 1.0);
//    assert(quiescentConcentration >= 0.0);
//    mQuiescentConcentration = quiescentConcentration;
//}
//
//double UniformCellCycleModelWithQuiescence::GetCriticalHypoxicDuration() const
//{
//    return mCriticalHypoxicDuration;
//}
//
//void UniformCellCycleModelWithQuiescence::SetCriticalHypoxicDuration(double criticalHypoxicDuration)
//{
//    assert(criticalHypoxicDuration >= 0.0);
//    mCriticalHypoxicDuration = criticalHypoxicDuration;
//}
//
//void UniformCellCycleModelWithQuiescence::SetCurrentHypoxiaOnsetTime(double currentHypoxiaOnsetTime)
//{
//    assert(currentHypoxiaOnsetTime >= 0.0);
//    mCurrentHypoxiaOnsetTime = currentHypoxiaOnsetTime;
//}
//
//// Overridden ReadyToDivide
//bool UniformCellCycleModelWithQuiescence::ReadyToDivide()
//{
//    assert(mpCell != nullptr);
//    UpdateHypoxicDuration();
//
//    if (!mReadyToDivide)
//    {
//        if (GetAge() >= mCellCycleDuration )
//        {
//            mReadyToDivide = true;
//        }
//    }
//    return mReadyToDivide;
//}
//
//
//void UniformCellCycleModelWithQuiescence::OutputCellCycleModelParameters(out_stream& rParamsFile)
//{
//    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
//    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";
//    *rParamsFile << "\t\t\t<HypoxicConcentration>" << mHypoxicConcentration << "</HypoxicConcentration>\n";
//    *rParamsFile << "\t\t\t<QuiescentConcentration>" << mQuiescentConcentration << "</QuiescentConcentration>\n";
//    *rParamsFile << "\t\t\t<CriticalHypoxicDuration>" << mCriticalHypoxicDuration << "</CriticalHypoxicDuration>\n";
//
//    // Call method on direct parent class
//    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
//}
//
//// Serialization for Boost >= 1.36
//#include "SerializationExportWrapperForCpp.hpp"
//CHASTE_CLASS_EXPORT(UniformCellCycleModelWithQuiescence)

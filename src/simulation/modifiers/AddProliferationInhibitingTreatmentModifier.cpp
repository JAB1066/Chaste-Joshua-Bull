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

#include "AddProliferationInhibitingTreatmentModifier.hpp"
#include "NearMacrophageCellProperty.hpp"
#include "UniformContactInhibitionCellCycleModelSlowByFibroblasts.hpp"

template<unsigned DIM>
AddProliferationInhibitingTreatmentModifier<DIM>::AddProliferationInhibitingTreatmentModifier()
: AbstractCellBasedSimulationModifier<DIM>(),
  mMinCellCycleDuration_PostTreatment(DBL_MAX),
  mMaxCellCycleDuration_PostTreatment(DBL_MAX),
  mTreatmentStartTime(DBL_MAX)
  {
  }

template<unsigned DIM>
AddProliferationInhibitingTreatmentModifier<DIM>::~AddProliferationInhibitingTreatmentModifier()
{
}

template<unsigned DIM>
double AddProliferationInhibitingTreatmentModifier<DIM>::GetTreatmentStartTime()
{
	return mTreatmentStartTime;
}

template<unsigned DIM>
void AddProliferationInhibitingTreatmentModifier<DIM>::SetTreatmentStartTime(double newTreatmentStartTime)
{
	mTreatmentStartTime = newTreatmentStartTime;
}

template<unsigned DIM>
double AddProliferationInhibitingTreatmentModifier<DIM>::GetMinCellCycleDuration_PostTreatment()
{
	return mMinCellCycleDuration_PostTreatment;
}

template<unsigned DIM>
void AddProliferationInhibitingTreatmentModifier<DIM>::SetMinCellCycleDuration_PostTreatment(double newMinCellCycleDuration_PostTreatment)
{
	mMinCellCycleDuration_PostTreatment = newMinCellCycleDuration_PostTreatment;
}

template<unsigned DIM>
double AddProliferationInhibitingTreatmentModifier<DIM>::GetMaxCellCycleDuration_PostTreatment()
{
	return mMaxCellCycleDuration_PostTreatment;
}

template<unsigned DIM>
void AddProliferationInhibitingTreatmentModifier<DIM>::SetMaxCellCycleDuration_PostTreatment(double newMaxCellCycleDuration_PostTreatment)
{
	mMaxCellCycleDuration_PostTreatment = newMaxCellCycleDuration_PostTreatment;
}

template<unsigned DIM>
double AddProliferationInhibitingTreatmentModifier<DIM>::GetMinCellCycleDuration_ByFibroblastPostTreatment()
{
	return mMinCellCycleDuration_ByFibroblastPostTreatment;
}

template<unsigned DIM>
void AddProliferationInhibitingTreatmentModifier<DIM>::SetMinCellCycleDuration_ByFibroblastPostTreatment(double newMinCellCycleDuration_ByFibroblastPostTreatment)
{
	mMinCellCycleDuration_ByFibroblastPostTreatment = newMinCellCycleDuration_ByFibroblastPostTreatment;
}

template<unsigned DIM>
double AddProliferationInhibitingTreatmentModifier<DIM>::GetMaxCellCycleDuration_ByFibroblastPostTreatment()
{
	return mMaxCellCycleDuration_ByFibroblastPostTreatment;
}

template<unsigned DIM>
void AddProliferationInhibitingTreatmentModifier<DIM>::SetMaxCellCycleDuration_ByFibroblastPostTreatment(double newMaxCellCycleDuration_ByFibroblastPostTreatment)
{
	mMaxCellCycleDuration_ByFibroblastPostTreatment = newMaxCellCycleDuration_ByFibroblastPostTreatment;
}

template<unsigned DIM>
void AddProliferationInhibitingTreatmentModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void AddProliferationInhibitingTreatmentModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void AddProliferationInhibitingTreatmentModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // For each cell in the cell population, we check if it's labelled as being close to a fibroblast ("macrophage").
    // We alter the proliferation rate accordingly.
    if (SimulationTime::Instance()->GetTime() == mTreatmentStartTime)
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
            // We access the cell cycle model for that cell and update the cell cycle
            double mutationState = cell_iter->GetMutationState()->GetColour();
            if(mutationState == 0) // If Sensitive
            {
                UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = static_cast<UniformContactInhibitionCellCycleModelSlowByFibroblasts*>(cell_iter->GetCellCycleModel());
                p_model->SetMinCellCycleDuration(mMinCellCycleDuration_PostTreatment);
                p_model->SetMaxCellCycleDuration(mMaxCellCycleDuration_PostTreatment);
                p_model->SetMinCellCycleDurationByFibroblast(mMinCellCycleDuration_ByFibroblastPostTreatment);
                p_model->SetMinCellCycleDurationByFibroblast(mMaxCellCycleDuration_ByFibroblastPostTreatment);
            }
            //UniformContactInhibitionCellCycleModelSlowByFibroblasts<2> cellCycleModel = cell_iter->GetCellCycleModel();
            //auto cellCycleModel = static_cast<ContactInhibitionCellCycleModel<2>* >(cell_iter->GetCellCycleModel());

            //cell_iter->GetCellCycleModel()->template SetTransitCellG1Duration(mProximalProliferationG1Duration);
        }
    }

	// Make sure the cell population is updated
	rCellPopulation.Update();

}

template<unsigned DIM>
void AddProliferationInhibitingTreatmentModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
	// No parameters to output, so just call method on direct parent class
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AddProliferationInhibitingTreatmentModifier<1>;
template class AddProliferationInhibitingTreatmentModifier<2>;
template class AddProliferationInhibitingTreatmentModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AddProliferationInhibitingTreatmentModifier)


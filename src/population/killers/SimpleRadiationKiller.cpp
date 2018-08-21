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

#include "SimpleRadiationKiller.hpp"
#include "SimulationTime.hpp"

template<unsigned DIM>
SimpleRadiationKiller<DIM>::SimpleRadiationKiller(AbstractCellPopulation<DIM>* pCellPopulation, double probabilityOfDeathInAnHour)
: AbstractCellKiller<DIM>(pCellPopulation),
  mProbabilityOfDeathInAnHour(probabilityOfDeathInAnHour)
  {
	if ((mProbabilityOfDeathInAnHour<0) || (mProbabilityOfDeathInAnHour>1))
	{
		EXCEPTION("Probability of death must be between zero and one");
	}
  }

template<unsigned DIM>
double SimpleRadiationKiller<DIM>::GetDeathProbabilityInAnHour() const
{
	return mProbabilityOfDeathInAnHour;
}

template<unsigned DIM>
void SimpleRadiationKiller<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
	/*
	 * We assume a constant time step and that there are an integer number (n = 1/dt)
	 * of time steps per hour. We also assume that this method is called every time step
	 * and that the probabilities of dying at different times are independent.
	 *
	 * Let q=mProbabilityOfDeathInAnHour and p="probability of death in a given time step".
	 *
	 * Probability of not dying in an hour:
	 * (1-q) = (1-p)^n = (1-p)^(1/dt).
	 *
	 * Rearranging for p:
	 * p = 1 - (1-q)^dt.
	 */
	double oxygenConc = pCell->GetCellData()->GetItem("oxygen");
	double death_prob_this_timestep = 1.0 - pow((1.0 - 2*oxygenConc*mProbabilityOfDeathInAnHour), SimulationTime::Instance()->GetTimeStep());
	// Oxygen enhances impact of radiation, so we scale the probability of death by the oxygen concentration
	if (!pCell->HasApoptosisBegun() &&
			RandomNumberGenerator::Instance()->ranf() < death_prob_this_timestep)
	{
		pCell->StartApoptosis();
	}
}

template<unsigned DIM>
void SimpleRadiationKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
	double time = SimulationTime::Instance()->GetTime();
	// Radiation times: 100,120,140,160,180,200 hours
	if (std::remainder(std::floor(time), 20) == 0 and time > 99 and time < 201)
	{
		for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
				cell_iter != this->mpCellPopulation->End();
				++cell_iter)
		{
			CheckAndLabelSingleCellForApoptosis(*cell_iter);
		}
	}
}

template<unsigned DIM>
void SimpleRadiationKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<ProbabilityOfDeathInAnHour>" << mProbabilityOfDeathInAnHour << "</ProbabilityOfDeathInAnHour>\n";

	// Call method on direct parent class
	AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class SimpleRadiationKiller<1>;
template class SimpleRadiationKiller<2>;
template class SimpleRadiationKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimpleRadiationKiller)

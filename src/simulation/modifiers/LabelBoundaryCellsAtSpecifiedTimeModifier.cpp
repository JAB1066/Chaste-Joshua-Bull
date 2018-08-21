/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "LabelBoundaryCellsAtSpecifiedTimeModifier.hpp"
#include "SmartPointers.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellLabel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Debug.hpp"
#include "RandomNumberGenerator.hpp"
#include "NodesOnlyMesh.hpp"

template<unsigned DIM>
LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::LabelBoundaryCellsAtSpecifiedTimeModifier()
: AbstractCellBasedSimulationModifier<DIM>(),
  timeToLabelCells(0.0),
  numberOfCellsToLabel(100),
  haveCellsBeenLabelled(false),
  numberOfHoursToRunSimulationAfterLabellingCells(UINT_MAX)
  {
  }

template<unsigned DIM>
LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::~LabelBoundaryCellsAtSpecifiedTimeModifier()
{
}

template<unsigned DIM>
void LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
double LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::GetTimeToLabelCells()
{
	return timeToLabelCells;
}

template<unsigned DIM>
void LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::SetTimeToLabelCells(double newTimeToLabelCells)
{
	timeToLabelCells = newTimeToLabelCells;
}

template<unsigned DIM>
unsigned LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::GetNumberOfCellsToLabel()
{
	return numberOfCellsToLabel;
}

template<unsigned DIM>
void LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::SetNumberOfCellsToLabel(unsigned newNumberOfCellsToLabel)
{
	numberOfCellsToLabel = newNumberOfCellsToLabel;
}

template<unsigned DIM>
unsigned LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::GetNumberOfHoursToRunSimulationAfterLabellingCells()
{
	return numberOfHoursToRunSimulationAfterLabellingCells;
}

template<unsigned DIM>
void LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::SetNumberOfHoursToRunSimulationAfterLabellingCells(unsigned newNumberOfHoursToRunSimulationAfterLabellingCells)
{
	numberOfHoursToRunSimulationAfterLabellingCells = newNumberOfHoursToRunSimulationAfterLabellingCells;
}

template<unsigned DIM>
void LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	UpdateCellData(rCellPopulation);
}



template<unsigned DIM>
void LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Make sure the cell population is updated
	rCellPopulation.Update();
	c_vector<double, DIM> centroid = rCellPopulation.GetCentroidOfCellPopulation();

	double time = SimulationTime::Instance()->GetTime();

	if((not haveCellsBeenLabelled) and time >= timeToLabelCells)
	{
		auto nbcp = dynamic_cast<NodeBasedCellPopulation<DIM>* >(&rCellPopulation);
		unsigned numBoundaryCells = 0;
		for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
				node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
				++node_iter)
		{
			if(node_iter->IsBoundaryNode())
			{
				numBoundaryCells++;
			}
		}

		// If numBoundaryCells < numberOfCellsToLabel, we label all of them...if not, we randomly select numberOfCellsToLabel from the list and label them
		MAKE_PTR(CellLabel, p_label);
		if(numBoundaryCells < numberOfCellsToLabel)
		{
			// Label all boundary cells
			for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
					node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
					++node_iter)
			{
				if(node_iter->IsBoundaryNode())
				{
					unsigned index = node_iter->GetIndex();
					CellPtr p_cell = nbcp->GetCellUsingLocationIndex(index);
					p_cell->AddCellProperty(p_label);
				}
			}
		}
		else
		{
			// We select numberOfCellsToLabel indices from 0 to numBoundaryCells
			std::vector<unsigned> shuffled_indices;
			RandomNumberGenerator::Instance()->Shuffle(numBoundaryCells, shuffled_indices);

			unsigned boundaryCellIndex = 0;
			unsigned numLabelledCells = 0;
			// Label the first numberOfCellsToLabel boundary cells
			for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
					node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
					++node_iter)
			{
				if(node_iter->IsBoundaryNode())
				{
					bool found = false;
					for(unsigned i = 0; i < numberOfCellsToLabel;i++)
					{
						if (shuffled_indices[i] == boundaryCellIndex)
						{
							// The i-th boundary cell has been chosen as one of the ones to be labelled
							found = true;
							break;
						}
					}

					if(found)
					{
						// Then we label the cell
						unsigned index = node_iter->GetIndex();
						CellPtr p_cell = nbcp->GetCellUsingLocationIndex(index);
						p_cell->AddCellProperty(p_label);
						numLabelledCells++;
					}
					// Increment the count of boundary cells!
					boundaryCellIndex++;
				}
				if(numLabelledCells > numberOfCellsToLabel)
				{
					break;
				}
			}
		}


		// Now ensure that we don't label cells again
		haveCellsBeenLabelled = true;

		// Set simulation ending time
		if(numberOfHoursToRunSimulationAfterLabellingCells < UINT_MAX)
		{
			double idealTimestep = SimulationTime::Instance()->GetTimeStep();

			SimulationTime::Destroy();

			double newEndTime = time + numberOfHoursToRunSimulationAfterLabellingCells;
			unsigned remainingTimesteps = numberOfHoursToRunSimulationAfterLabellingCells/idealTimestep+0.5;
			unsigned totalTimesteps = newEndTime/idealTimestep+0.5;
			PRINT_VARIABLE(time);
			PRINT_VARIABLE(newEndTime);
			PRINT_VARIABLE(remainingTimesteps);


			SimulationTime::Instance()->SetStartTime(time);
			SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(newEndTime,remainingTimesteps);
			SimulationTime::Instance()->SetStartTime(time);
		}
	}

}

template<unsigned DIM>
void LabelBoundaryCellsAtSpecifiedTimeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

	*rParamsFile << "\t\t\t<timeToLabelCells>" << timeToLabelCells << "</timeToLabelCells>\n";
	*rParamsFile << "\t\t\t<numberOfCellsToLabel>" << numberOfCellsToLabel << "</numberOfCellsToLabel>\n";
	*rParamsFile << "\t\t\t<numberOfHoursToRunSimulationAfterLabellingCells>" << numberOfHoursToRunSimulationAfterLabellingCells << "</numberOfHoursToRunSimulationAfterLabellingCells>\n";


	// Next, call method on direct parent class
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class LabelBoundaryCellsAtSpecifiedTimeModifier<1>;
template class LabelBoundaryCellsAtSpecifiedTimeModifier<2>;
template class LabelBoundaryCellsAtSpecifiedTimeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LabelBoundaryCellsAtSpecifiedTimeModifier)


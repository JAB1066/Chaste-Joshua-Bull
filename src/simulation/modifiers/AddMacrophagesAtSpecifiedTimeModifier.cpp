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

#include "AddMacrophagesAtSpecifiedTimeModifier.hpp"
#include "SmartPointers.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Debug.hpp"
#include "RandomNumberGenerator.hpp"
#include "NodesOnlyMesh.hpp"

template<unsigned DIM>
AddMacrophagesAtSpecifiedTimeModifier<DIM>::AddMacrophagesAtSpecifiedTimeModifier()
: AbstractCellBasedSimulationModifier<DIM>(),
  timeToAddMacrophages(0.0),
  numberOfMacrophagesToAdd(100),
  haveMacrophagesBeenAdded(false),
  numberOfHoursToRunSimulationAfterAddingMacrophages(UINT_MAX)
  {
  }

template<unsigned DIM>
AddMacrophagesAtSpecifiedTimeModifier<DIM>::~AddMacrophagesAtSpecifiedTimeModifier()
{
}

template<unsigned DIM>
void AddMacrophagesAtSpecifiedTimeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
double AddMacrophagesAtSpecifiedTimeModifier<DIM>::GetTimeToAddMacrophages()
{
	return timeToAddMacrophages;
}

template<unsigned DIM>
void AddMacrophagesAtSpecifiedTimeModifier<DIM>::SetTimeToAddMacrophages(double newTimeToAddMacrophages)
{
	timeToAddMacrophages = newTimeToAddMacrophages;
}

template<unsigned DIM>
unsigned AddMacrophagesAtSpecifiedTimeModifier<DIM>::GetNumberOfMacrophagesToAdd()
{
	return numberOfMacrophagesToAdd;
}

template<unsigned DIM>
void AddMacrophagesAtSpecifiedTimeModifier<DIM>::SetNumberOfMacrophagesToAdd(unsigned newNumberOfMacrophagesToAdd)
{
	numberOfMacrophagesToAdd = newNumberOfMacrophagesToAdd;
}

template<unsigned DIM>
unsigned AddMacrophagesAtSpecifiedTimeModifier<DIM>::GetNumberOfHoursToRunSimulationAfterAddingMacrophages()
{
	return numberOfHoursToRunSimulationAfterAddingMacrophages;
}

template<unsigned DIM>
void AddMacrophagesAtSpecifiedTimeModifier<DIM>::SetNumberOfHoursToRunSimulationAfterAddingMacrophages(unsigned newNumberOfHoursToRunSimulationAfterAddingMacrophages)
{
	numberOfHoursToRunSimulationAfterAddingMacrophages = newNumberOfHoursToRunSimulationAfterAddingMacrophages;
}

template<unsigned DIM>
void AddMacrophagesAtSpecifiedTimeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	UpdateCellData(rCellPopulation);
}



template<unsigned DIM>
void AddMacrophagesAtSpecifiedTimeModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Make sure the cell population is updated
	rCellPopulation.Update();

	double time = SimulationTime::Instance()->GetTime();

	if((not haveMacrophagesBeenAdded) and time >= timeToAddMacrophages)
	{
		auto nbcp = dynamic_cast<NodeBasedCellPopulation<DIM>* >(&rCellPopulation);

		/*
		 * Now we add our macrophages
		 */
		//unsigned nodeNum = rCellPopulation.GetNumNodes();

		std::list<CellPtr> tumour_cells = rCellPopulation.rGetCells();
		std::list<CellPtr>::iterator it;
		double macrophageSphereRadius = 0;
		for (it = tumour_cells.begin(); it != tumour_cells.end(); it++)
		{
			double radius = 0;
			c_vector<double, DIM> location = rCellPopulation.GetLocationOfCellCentre(*it);
			for (unsigned i=0; i<DIM; i++)
			{
				radius += pow(location[i],2);
			}
			if(radius > macrophageSphereRadius)
			{
				macrophageSphereRadius = radius;
			}
		}

		// Macrophage Nodes - in a shell at fixed radius
		// Number of Macrophages remains constant throughout simulation

		unsigned macCount=0;
		double random_x;
		double random_y;
		double random_z;
		double normalize;
		macrophageSphereRadius =  std::sqrt(macrophageSphereRadius) + 0.5;


		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);
		MAKE_PTR(WildTypeCellMutationState, p_state);

		// Add nodes at random points at edge of spheroid
		while(macCount < numberOfMacrophagesToAdd)
		{

			random_x = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
			random_y = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
			if(DIM == 2)
			{
				normalize = macrophageSphereRadius/sqrt(pow(random_x,2) + pow(random_y,2));
			}
			if(DIM == 3)
			{
				random_z = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
				normalize = macrophageSphereRadius/sqrt(pow(random_x,2) + pow(random_y,2)+ pow(random_z,2));
			}

			// Make vectors lie on sphere of radiusmacrophageSphereRadius

			// Make Macrophage
			NoCellCycleModel* p_model = new NoCellCycleModel;
			p_model->SetDimension(DIM);
			CellPtr pNewCell(new Cell(p_state, p_model));
			pNewCell->SetCellProliferativeType(p_macrophage_type);
			pNewCell->GetCellData()->SetItem("oxygen", 1);

			//node = new Node<3>(nodeNum,  false,  random_x*normalize, random_y*normalize, random_z*normalize);

			// 10000 here is a placeholder - it should be overwritten when we add the node to the nbcp...I hope
			Node<DIM>* p_new_node;
			if(DIM == 1)
			{
				p_new_node = new Node<DIM>(100000,  false,  random_x*normalize);
			}
			if(DIM == 2)
			{
				p_new_node = new Node<DIM>(100000,  false,  random_x*normalize, random_y*normalize);
			}
			if(DIM == 3)
			{
				p_new_node = new Node<DIM>(100000,  false,  random_x*normalize, random_y*normalize, random_z*normalize);
			}

			p_new_node->ClearAppliedForce(); // Incase velocity is ouptut on the same timestep as the cell has divided
			//p_new_node->SetRadius(0.05); // Beads are much smaller than cells
			//NodesOnlyMesh<3> mesh = nbcp->rGetMesh();
			unsigned new_node_index = nbcp->rGetMesh().AddNode(p_new_node);


			// Update cells vector
			nbcp->rGetCells().push_back(pNewCell);

			// Update mappings between cells and location indices
			nbcp->SetCellUsingLocationIndex(new_node_index, pNewCell);

			macCount++;
		}

		// Now ensure that we don't add macrophages again
		haveMacrophagesBeenAdded = true;

		// Set simulation ending time
		if(numberOfHoursToRunSimulationAfterAddingMacrophages < UINT_MAX)
		{
			double idealTimestep = SimulationTime::Instance()->GetTimeStep();

			SimulationTime::Destroy();

			double newEndTime = time + numberOfHoursToRunSimulationAfterAddingMacrophages;
			unsigned remainingTimesteps = numberOfHoursToRunSimulationAfterAddingMacrophages/idealTimestep+0.5;
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
void AddMacrophagesAtSpecifiedTimeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

	*rParamsFile << "\t\t\t<timeToAddMacrophages>" << timeToAddMacrophages << "</timeToAddMacrophages>\n";
	*rParamsFile << "\t\t\t<numberOfMacrophagesToAdd>" << numberOfMacrophagesToAdd << "</numberOfMacrophagesToAdd>\n";
	*rParamsFile << "\t\t\t<numberOfHoursToRunSimulationAfterAddingMacrophages>" << numberOfHoursToRunSimulationAfterAddingMacrophages << "</numberOfHoursToRunSimulationAfterAddingMacrophages>\n";


	// Next, call method on direct parent class
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AddMacrophagesAtSpecifiedTimeModifier<1>;
template class AddMacrophagesAtSpecifiedTimeModifier<2>;
template class AddMacrophagesAtSpecifiedTimeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AddMacrophagesAtSpecifiedTimeModifier)


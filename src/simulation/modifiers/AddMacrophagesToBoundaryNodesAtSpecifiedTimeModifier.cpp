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

#include "AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier.hpp"
#include "SmartPointers.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Debug.hpp"
#include "RandomNumberGenerator.hpp"
#include "NodesOnlyMesh.hpp"

template<unsigned DIM>
AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier()
: AbstractCellBasedSimulationModifier<DIM>(),
  timeToAddMacrophages(0.0),
  numberOfMacrophagesToAdd(100),
  haveMacrophagesBeenAdded(false),
  numberOfHoursToRunSimulationAfterAddingMacrophages(UINT_MAX)
  {
  }

template<unsigned DIM>
AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::~AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier()
{
}

template<unsigned DIM>
void AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
double AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::GetTimeToAddMacrophages()
{
	return timeToAddMacrophages;
}

template<unsigned DIM>
void AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::SetTimeToAddMacrophages(double newTimeToAddMacrophages)
{
	timeToAddMacrophages = newTimeToAddMacrophages;
}

template<unsigned DIM>
unsigned AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::GetNumberOfMacrophagesToAdd()
{
	return numberOfMacrophagesToAdd;
}

template<unsigned DIM>
void AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::SetNumberOfMacrophagesToAdd(unsigned newNumberOfMacrophagesToAdd)
{
	numberOfMacrophagesToAdd = newNumberOfMacrophagesToAdd;
}

template<unsigned DIM>
unsigned AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::GetNumberOfHoursToRunSimulationAfterAddingMacrophages()
{
	return numberOfHoursToRunSimulationAfterAddingMacrophages;
}

template<unsigned DIM>
void AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::SetNumberOfHoursToRunSimulationAfterAddingMacrophages(unsigned newNumberOfHoursToRunSimulationAfterAddingMacrophages)
{
	numberOfHoursToRunSimulationAfterAddingMacrophages = newNumberOfHoursToRunSimulationAfterAddingMacrophages;
}

template<unsigned DIM>
void AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	UpdateCellData(rCellPopulation);
}



template<unsigned DIM>
void AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Make sure the cell population is updated
	rCellPopulation.Update();
	c_vector<double, DIM> centroid = rCellPopulation.GetCentroidOfCellPopulation();

	double time = SimulationTime::Instance()->GetTime();

	if((not haveMacrophagesBeenAdded) and time >= timeToAddMacrophages)
	{
		auto nbcp = dynamic_cast<NodeBasedCellPopulation<DIM>* >(&rCellPopulation);

		std::vector<c_vector<double, DIM> > boundaryCellLocations;
		for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
				node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
				++node_iter)
		{
			if(node_iter->IsBoundaryNode())
			{
				c_vector<double, DIM> location = node_iter->rGetLocation();
				boundaryCellLocations.push_back(location);
			}
		}


		/*
		 * We randomly select boundary nodes, and try to add a macrophage somewhere around it on the outside of the tumour
		 */
		unsigned macCount=0;
		double random_x;
		double random_y;
		double random_z;
		double normalize;
		while(macCount < numberOfMacrophagesToAdd)
		{
			// Choose a random boundary location
			c_vector<double, DIM> boundaryLocation = boundaryCellLocations[RandomNumberGenerator::Instance()->randMod(boundaryCellLocations.size())];

			// Make a perturbation
			c_vector<double, DIM> perturbation;

			random_x = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
			random_y = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
			if(DIM == 2)
			{
				normalize = 1/sqrt(pow(random_x,2) + pow(random_y,2));
			}
			if(DIM == 3)
			{
				random_z = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
				normalize = 1/sqrt(pow(random_x,2) + pow(random_y,2)+ pow(random_z,2));
			}

			// Calculate vector away from centre (n)
			double perturbationSize = 0;
			double normalSize = 0;
			double normalDotPerturbation = 0;
			c_vector<double, DIM> normal;
			for (unsigned i=0; i<DIM; i++)
			{
				perturbation[i] = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
				normal[i] = boundaryLocation[i] - centroid[i];

				normalDotPerturbation += perturbation[i]*normal[i];

				normalSize += pow(normal[i],2);
				perturbationSize += pow(perturbation[i],2);
			}
			perturbationSize = sqrt(perturbationSize);
			normalSize = sqrt(normalSize);

			// If the angle between the perturbation and the normal is less than pi/2, then we can try to place the macrophage (as outside the sphere)
			double theta = normalDotPerturbation/(perturbationSize*normalSize);
			if(theta < M_PI/2)
			{
				MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);
				MAKE_PTR(WildTypeCellMutationState, p_state);

				// Make Macrophage
				NoCellCycleModel* p_model = new NoCellCycleModel;
				p_model->SetDimension(DIM);
				CellPtr pNewCell(new Cell(p_state, p_model));
				pNewCell->SetCellProliferativeType(p_macrophage_type);
				pNewCell->GetCellData()->SetItem("oxygen", 1);


				// 10000 here is a placeholder - it should be overwritten when we add the node to the nbcp...I hope
				Node<DIM>* p_new_node;
				if(DIM == 1)
				{
					p_new_node = new Node<DIM>(100000,  false,  boundaryLocation[0]+perturbation[0]/perturbationSize);
				}
				if(DIM == 2)
				{
					p_new_node = new Node<DIM>(100000,  false,  boundaryLocation[0]+perturbation[0]/perturbationSize, boundaryLocation[1]+perturbation[1]/perturbationSize);
				}
				if(DIM == 3)
				{
					p_new_node = new Node<DIM>(100000,  false,  boundaryLocation[0]+perturbation[0]/perturbationSize, boundaryLocation[1]+perturbation[1]/perturbationSize, boundaryLocation[2]+perturbation[2]/perturbationSize);
				}

				p_new_node->ClearAppliedForce(); // Incase velocity is ouptut on the same timestep as the cell has divided
				p_new_node->SetRadius(0.5*0.75); // Beads are much smaller than cells
				//NodesOnlyMesh<3> mesh = nbcp->rGetMesh();
				unsigned new_node_index = nbcp->rGetMesh().AddNode(p_new_node);

				// Update cells vector
				nbcp->rGetCells().push_back(pNewCell);

				// Update mappings between cells and location indices
				nbcp->SetCellUsingLocationIndex(new_node_index, pNewCell);

				macCount++;

			}

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
void AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

	*rParamsFile << "\t\t\t<timeToAddMacrophages>" << timeToAddMacrophages << "</timeToAddMacrophages>\n";
	*rParamsFile << "\t\t\t<numberOfMacrophagesToAdd>" << numberOfMacrophagesToAdd << "</numberOfMacrophagesToAdd>\n";
	*rParamsFile << "\t\t\t<numberOfHoursToRunSimulationAfterAddingMacrophages>" << numberOfHoursToRunSimulationAfterAddingMacrophages << "</numberOfHoursToRunSimulationAfterAddingMacrophages>\n";


	// Next, call method on direct parent class
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<1>;
template class AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<2>;
template class AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier)


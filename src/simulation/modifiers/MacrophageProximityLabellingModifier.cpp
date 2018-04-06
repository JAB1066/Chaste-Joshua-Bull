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

#include "MacrophageProximityLabellingModifier.hpp"
#include "SmartPointers.hpp"
#include "NearMacrophageCellProperty.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned DIM>
MacrophageProximityLabellingModifier<DIM>::MacrophageProximityLabellingModifier()
: AbstractCellBasedSimulationModifier<DIM>(),
  macrophageProximityLabellingRadius(1.0)
  {
  }

template<unsigned DIM>
MacrophageProximityLabellingModifier<DIM>::~MacrophageProximityLabellingModifier()
{
}

template<unsigned DIM>
void MacrophageProximityLabellingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
double MacrophageProximityLabellingModifier<DIM>::GetMacrophageProximityLabellingRadius()
{
	return macrophageProximityLabellingRadius;
}

template<unsigned DIM>
void MacrophageProximityLabellingModifier<DIM>::SetMacrophageProximityLabellingRadius(double newMacrophageProximityLabellingRadius)
{
	assert(newMacrophageProximityLabellingRadius > 0);
	macrophageProximityLabellingRadius = newMacrophageProximityLabellingRadius;
}

template<unsigned DIM>
void MacrophageProximityLabellingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
double MacrophageProximityLabellingModifier<DIM>::CalculateDistanceFromMacrophage(unsigned nodeGlobalIndex, AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Get node location
	Node<DIM>* p_node = rCellPopulation.GetNode(nodeGlobalIndex);
	const c_vector<double, DIM>& r_node_location = p_node->rGetLocation();
	//Node<SPACE_DIM>* p_node = rCellPopulation.GetNode(nodeGlobalIndex);
	//const c_vector<double, SPACE_DIM>& r_node_location = p_node->rGetLocation();

	//std::vector<double> distancesVector;
	double minDist = DBL_MAX;

	// Iterate over cell population to find macrophages
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
			cell_iter != rCellPopulation.End();
			++cell_iter)
	{
	// If cell is a macrophage...
		//...calculate distance between cells
		if (cell_iter->GetCellProliferativeType()->template IsType<MacrophageCellProliferativeType>())
		{
			// Get node location
			unsigned node_b_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
			Node<DIM>* p_node_b = rCellPopulation.GetNode(node_b_index);


			const c_vector<double, DIM>& r_node_b_location = p_node_b->rGetLocation();

			// Get vector separating nodes
			c_vector<double, DIM> difference;
			/*
			 * We use the mesh method GetVectorFromAtoB() to compute the direction of the
			 * unit vector along the line joining the two nodes, rather than simply subtract
			 * their positions, because this method can be overloaded (e.g. to enforce a
			 * periodic boundary in Cylindrical2dMesh).
			 */
			difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_location, r_node_b_location);

			if (minDist > norm_2(difference))
            {
                minDist = norm_2(difference);
            }

			// Calculate the distance between the two nodes
			//double distance_between_nodes = norm_2(difference);
			//assert(!std::isnan(distance_between_nodes));

			//distancesVector.push_back(distance_between_nodes);
		}
	}
	return minDist;//*std::min_element(distancesVector.begin(), distancesVector.end()); // Return smallest distance in vector
}


template<unsigned DIM>
void MacrophageProximityLabellingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Make sure the cell population is updated
	rCellPopulation.Update();

	MAKE_PTR(NearMacrophageCellProperty, p_MacrophageProximityLabel);

	// Iterate over cell population
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
			cell_iter != rCellPopulation.End();
			++cell_iter)
	{
		double distanceToMacrophage = CalculateDistanceFromMacrophage(rCellPopulation.GetLocationIndexUsingCell(*cell_iter), rCellPopulation);

		// If close to macrophage, label the cell
		if(distanceToMacrophage < macrophageProximityLabellingRadius)
		{
			if(!(cell_iter->template HasCellProperty<NearMacrophageCellProperty>()))
			{
				cell_iter->AddCellProperty(p_MacrophageProximityLabel);;
			}
		}
		else
		{
			if(cell_iter->template HasCellProperty<NearMacrophageCellProperty>())
			{
				cell_iter->template RemoveCellProperty<NearMacrophageCellProperty>();
			}
		}
	}
}

template<unsigned DIM>
void MacrophageProximityLabellingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

    *rParamsFile << "\t\t\t<MacrophageProximityLabellingRadius>" << macrophageProximityLabellingRadius << "</MacrophageProximityLabellingRadius>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MacrophageProximityLabellingModifier<1>;
template class MacrophageProximityLabellingModifier<2>;
template class MacrophageProximityLabellingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MacrophageProximityLabellingModifier)


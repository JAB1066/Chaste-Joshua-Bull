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

#include "VaryDampingCoefficientOfBoundaryNodes.hpp"

template<unsigned DIM>
VaryDampingCoefficientOfBoundaryNodes<DIM>::VaryDampingCoefficientOfBoundaryNodes()
: AbstractForce<DIM>(),
  mBoundaryNodeDampingConstant(1.0),
  mInteriorNodeDampingConstant(1.0)
  {
  }

template<unsigned DIM>
VaryDampingCoefficientOfBoundaryNodes<DIM>::~VaryDampingCoefficientOfBoundaryNodes()
{
}

template<unsigned DIM>
double VaryDampingCoefficientOfBoundaryNodes<DIM>::GetBoundaryNodeDampingConstant()
{
	return mBoundaryNodeDampingConstant;
}

template<unsigned DIM>
void VaryDampingCoefficientOfBoundaryNodes<DIM>::SetBoundaryNodeDampingConstant(double newBoundaryNodeDampingConstant)
{
	assert(newBoundaryNodeDampingConstant > 0);
	mBoundaryNodeDampingConstant = newBoundaryNodeDampingConstant;
}

template<unsigned DIM>
double VaryDampingCoefficientOfBoundaryNodes<DIM>::GetInteriorNodeDampingConstant()
{
	return mInteriorNodeDampingConstant;
}

template<unsigned DIM>
void VaryDampingCoefficientOfBoundaryNodes<DIM>::SetInteriorNodeDampingConstant(double newInteriorNodeDampingConstant)
{
	assert(newInteriorNodeDampingConstant > 0);
	mInteriorNodeDampingConstant = newInteriorNodeDampingConstant;
}

template<unsigned DIM>
void VaryDampingCoefficientOfBoundaryNodes<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

	c_vector<double, DIM> force_contribution;
	double dt = SimulationTime::Instance()->GetTimeStep();
	for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
			node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
			++node_iter)
	{

		// We make a hack to determine cells near the boundary. Since boundary nodes have oxygen concentration = 1.0,
		// any node which has that must be near the boundary. Yup, it's crude...can't find a better way, sorry!
		// This method is basically taken from AbstractBoxDomainPdeModifier

		// Get cell
		unsigned index = node_iter->GetIndex();
		CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(index);

		double oxygenLevel = p_cell->GetCellData()->GetItem("oxygen");

		// We scale only the component of the vector pointing OUT of the spheroid by 100 - any movement inwards remains unscaled.
		// To do this, for any boundary node we calculate whether or not the new force contribution moves it across the tangent plane at its point
		// If it does, scale. If it doesn't, OK.

		if( oxygenLevel > 0.99)
		{
			// Get the current applied force
			force_contribution = node_iter->rGetAppliedForce();

			c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(p_cell);
			if (inner_prod(cell_location + dt*force_contribution, cell_location) > 0.0)
			{

				// As node is on the boundary, scale its force by mBoundaryNodeDampingConstant.

				for (unsigned i=0; i<DIM; i++)
				{
					force_contribution[i] /= mBoundaryNodeDampingConstant;
				}

				// Now clear the current force and add the rescaled force contribution
				node_iter->ClearAppliedForce();
				node_iter->AddAppliedForceContribution(force_contribution);
			}
		}
		else
		{
			// Do nothing
			//node_iter->ClearAppliedForce();
		}





		//		// Working version
		//		if( oxygenLevel > 0.99)
		//		{
		//			// Get the current applied force
		//			force_contribution = node_iter->rGetAppliedForce();
		//
		//			// As node is on the boundary, scale its force by mBoundaryNodeDampingConstant.
		//
		//			for (unsigned i=0; i<DIM; i++)
		//			{
		//				force_contribution[i] /= mBoundaryNodeDampingConstant;
		//			}
		//
		//			// Now clear the current force and add the rescaled force contribution
		//			node_iter->ClearAppliedForce();
		//			node_iter->AddAppliedForceContribution(force_contribution);
		//		}
		//		else
		//		{
		//			// Do nothing
		//			//node_iter->ClearAppliedForce();
		//		}
	}

}

template<unsigned DIM>
void VaryDampingCoefficientOfBoundaryNodes<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<mBoundaryNodeDampingConstant>" << mBoundaryNodeDampingConstant << "</mBoundaryNodeDampingConstant> \n";
	*rParamsFile << "\t\t\t<mInteriorNodeDampingConstant>" << mInteriorNodeDampingConstant << "</mInteriorNodeDampingConstant> \n";

	// Call direct parent class
	AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class VaryDampingCoefficientOfBoundaryNodes<1>;
template class VaryDampingCoefficientOfBoundaryNodes<2>;
template class VaryDampingCoefficientOfBoundaryNodes<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VaryDampingCoefficientOfBoundaryNodes)

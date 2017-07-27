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

#include "CuboidPeriodicBoundaryCondition.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"
#include "SmartPointers.hpp"

template<unsigned DIM>
CuboidPeriodicBoundaryCondition<DIM>::CuboidPeriodicBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                                      ChastePoint<DIM> lower,
																	  ChastePoint<DIM> upper)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
      mLower(lower),
	  mUpper(upper)
{

	if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation) == NULL)
    {
        EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
    }
    if (DIM == 1)
    {
        EXCEPTION("This boundary condition is not implemented in 1D.");
    }
}

template<unsigned DIM>
ChastePoint<DIM> CuboidPeriodicBoundaryCondition<DIM>::GetLowerPoint() const
{
    return mLower;
}

template<unsigned DIM>
ChastePoint<DIM> CuboidPeriodicBoundaryCondition<DIM>::GetUpperPoint() const
{
    return mUpper;
}

template<unsigned DIM>
void CuboidPeriodicBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
    	bool moveCell = false;
        // If new location is outside cuboid in some dimension, move it to the other side
        c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        c_vector<double, DIM> new_cell_location = cell_location;
        for(unsigned i = 0; i<DIM; i++)
        {
        	if(cell_location[i] > mUpper[i])
        	{
        		new_cell_location[i] = mLower[i] + cell_location[i] - mUpper[i];
        		moveCell = true;
        	}
        	if(cell_location[i] < mLower[i])
        	{
        		new_cell_location[i] = mUpper[i] - (mLower[i] - cell_location[i]);
        		moveCell = true;
        	}
        }

        // If cell has moved, move it
        if(moveCell)
        {
        	unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        	Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

        	p_node->rGetModifiableLocation() = new_cell_location;
        }
    }
}

template<unsigned DIM>
bool CuboidPeriodicBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (mLower, mUpper));

    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
    	ChastePoint<DIM> cell_location(this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter));

    	if(not p_cuboid->DoesContain(cell_location))
    	{
            // ...then the boundary condition is not satisfied
            condition_satisfied = false;
            break;
    	}
    }
    return condition_satisfied;
}

template<unsigned DIM>
void CuboidPeriodicBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class CuboidPeriodicBoundaryCondition<1>;
template class CuboidPeriodicBoundaryCondition<2>;
template class CuboidPeriodicBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CuboidPeriodicBoundaryCondition)

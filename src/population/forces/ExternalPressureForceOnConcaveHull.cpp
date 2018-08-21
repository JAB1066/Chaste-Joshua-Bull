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

#include "ExternalPressureForceOnConcaveHull.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "AbstractCellProliferativeType.hpp"

template<unsigned DIM>
ExternalPressureForceOnConcaveHull<DIM>::ExternalPressureForceOnConcaveHull()
: AbstractForce<DIM>(),
  mPressure(1.0),
  mIsForceDirectedTowardsOrigin(true),
  mAlpha(1.0)
  {
  }

template<unsigned DIM>
ExternalPressureForceOnConcaveHull<DIM>::~ExternalPressureForceOnConcaveHull()
{
}

template<unsigned DIM>
double ExternalPressureForceOnConcaveHull<DIM>::GetPressure()
{
	return mPressure;
}

template<unsigned DIM>
void ExternalPressureForceOnConcaveHull<DIM>::SetPressure(double newPressure)
{
	assert(newPressure > 0);
	mPressure = newPressure;
}

template<unsigned DIM>
bool ExternalPressureForceOnConcaveHull<DIM>::GetIsForceDirectedTowardsOrigin()
{
	return mIsForceDirectedTowardsOrigin;
}

template<unsigned DIM>
void ExternalPressureForceOnConcaveHull<DIM>::SetIsForceDirectedTowardsOrigin(bool newIsForceDirectedTowardsOrigin)
{
	mIsForceDirectedTowardsOrigin = newIsForceDirectedTowardsOrigin;
}

template<unsigned DIM>
double ExternalPressureForceOnConcaveHull<DIM>::GetAlpha()
{
	return mAlpha;
}

template<unsigned DIM>
void ExternalPressureForceOnConcaveHull<DIM>::SetAlpha(double newAlpha)
{
	mAlpha = newAlpha;
}


template<unsigned DIM>
void ExternalPressureForceOnConcaveHull<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

	c_vector<double, DIM> force_contribution;

	c_vector<double, DIM> centroid;
	if(mIsForceDirectedTowardsOrigin) // Direct force from external pressure towards coordinate origin
	{
		centroid = zero_vector<double>(DIM);
	}
	else // Direct force from external pressure towards spheroid centroid
	{
		centroid = rCellPopulation.GetCentroidOfCellPopulation();
	}

	for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
			node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
			++node_iter)
	{
		// If node is a boundary node, add force in direction of centroid
		// Note - this is set in EllipticBoxDomainPdeModifierVariableTimestep - if you're not using that to find the alpha shape, it won't work!
		if( node_iter->IsBoundaryNode())
		{
			c_vector<double, DIM> cell_location = node_iter->rGetLocation();
			c_vector<double, DIM> vectorToCentroid = rCellPopulation.rGetMesh().GetVectorFromAtoB(cell_location, centroid);

			double cellRadialDistanceFromCentroid = 0;
			for (unsigned i=0; i<DIM; i++)
			{
				cellRadialDistanceFromCentroid+= vectorToCentroid[i]*vectorToCentroid[i];
			}
			cellRadialDistanceFromCentroid = sqrt(cellRadialDistanceFromCentroid);

			for (unsigned i=0; i<DIM; i++)
			{
				// We need force to be proportional to radius^2 to maintain density of cells
				//force_contribution[i] = cellRadialDistanceFromCentroid*cellRadialDistanceFromCentroid*( vectorToCentroid[i]/cellRadialDistanceFromCentroid ) * ( mPressure/pow(2*(node_iter->GetRadius()),2) ); // force = pressure / (cell diameter)^2
				force_contribution[i] = ( vectorToCentroid[i]/cellRadialDistanceFromCentroid ) * ( mPressure/pow(2*(node_iter->GetRadius()),2) ); // force = pressure / (cell diameter)^2
			}

			node_iter->AddAppliedForceContribution(force_contribution);

		}
	}

}

template<unsigned DIM>
void ExternalPressureForceOnConcaveHull<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<mPressure>" << mPressure << "</mPressure> \n";
	*rParamsFile << "\t\t\t<mAlpha>" << mAlpha << "</mAlpha> \n";

	// Call direct parent class
	AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ExternalPressureForceOnConcaveHull<1>;
template class ExternalPressureForceOnConcaveHull<2>;
template class ExternalPressureForceOnConcaveHull<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExternalPressureForceOnConcaveHull)

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

#include "ExternalPressureForceOnSurroundingSphere.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "AbstractCellProliferativeType.hpp"

template<unsigned DIM>
ExternalPressureForceOnSurroundingSphere<DIM>::ExternalPressureForceOnSurroundingSphere()
: AbstractForce<DIM>(),
  mPressure(1.0),
  mIsForceDirectedTowardsOrigin(true)
  {
  }

template<unsigned DIM>
ExternalPressureForceOnSurroundingSphere<DIM>::~ExternalPressureForceOnSurroundingSphere()
{
}

template<unsigned DIM>
double ExternalPressureForceOnSurroundingSphere<DIM>::GetPressure()
{
	return mPressure;
}

template<unsigned DIM>
void ExternalPressureForceOnSurroundingSphere<DIM>::SetPressure(double newPressure)
{
	assert(newPressure > 0);
	mPressure = newPressure;
}

template<unsigned DIM>
bool ExternalPressureForceOnSurroundingSphere<DIM>::GetIsForceDirectedTowardsOrigin()
{
	return mIsForceDirectedTowardsOrigin;
}

template<unsigned DIM>
void ExternalPressureForceOnSurroundingSphere<DIM>::SetIsForceDirectedTowardsOrigin(bool newIsForceDirectedTowardsOrigin)
{
	mIsForceDirectedTowardsOrigin = newIsForceDirectedTowardsOrigin;
}


template<unsigned DIM>
void ExternalPressureForceOnSurroundingSphere<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
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

	// Find max radius - this is not remotely efficient
	double currentMaxRadiusSquared = 0;
	for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
			node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
			++node_iter)
	{
		// Get cell
		unsigned index = node_iter->GetIndex();

		CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(index);
		boost::shared_ptr<AbstractCellProliferativeType> celltype = p_cell->GetCellProliferativeType();
		if (not celltype->template IsType<MacrophageCellProliferativeType>())
		{



			c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(p_cell);
			c_vector<double, DIM> force_contribution;

			c_vector<double, DIM> vectorToCentroid = rCellPopulation.rGetMesh().GetVectorFromAtoB(cell_location, centroid);

			double cellRadialDistanceFromCentroid = 0;
			for (unsigned i=0; i<DIM; i++)
			{
				cellRadialDistanceFromCentroid+= vectorToCentroid[i]*vectorToCentroid[i];
			}
			if (currentMaxRadiusSquared < cellRadialDistanceFromCentroid)
			{
				currentMaxRadiusSquared = cellRadialDistanceFromCentroid;
			}
		}
	}


	double currentMaxRadius = sqrt(currentMaxRadiusSquared);
	for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
			node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
			++node_iter)
	{

		// Look, I said this was inefficient, OK?! I'm just trying something out.

		// Get cell
		unsigned index = node_iter->GetIndex();

		CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(index);
		boost::shared_ptr<AbstractCellProliferativeType> celltype = p_cell->GetCellProliferativeType();
		//if (not celltype->template IsType<MacrophageCellProliferativeType>())
		//{

		// Pressure is force per (cell diameter)^2
		// Force is to be applied towards the origin, radially
		c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(p_cell);
		c_vector<double, DIM> vectorToCentroid = rCellPopulation.rGetMesh().GetVectorFromAtoB(cell_location, centroid);

		double cellRadialDistanceFromCentroid = 0;
		for (unsigned i=0; i<DIM; i++)
		{
			cellRadialDistanceFromCentroid+= vectorToCentroid[i]*vectorToCentroid[i];
		}
		cellRadialDistanceFromCentroid = sqrt(cellRadialDistanceFromCentroid);

		if( abs(currentMaxRadius - cellRadialDistanceFromCentroid) < 0.5) // If within 0.5 cell widths of "boundary"
		{
			c_vector<double, DIM> force_contribution;
			for (unsigned i=0; i<DIM; i++)
			{
				// We need force to be proportional to radius^2 to maintain density of cells
				//force_contribution[i] = cellRadialDistanceFromCentroid*cellRadialDistanceFromCentroid*( vectorToCentroid[i]/cellRadialDistanceFromCentroid ) * ( mPressure/pow(2*(node_iter->GetRadius()),2) ); // force = pressure / (cell diameter)^2
				force_contribution[i] = ( vectorToCentroid[i]/cellRadialDistanceFromCentroid ) * ( mPressure/pow(2*(node_iter->GetRadius()),2) ); // force = pressure / (cell diameter)^2
			}

			node_iter->AddAppliedForceContribution(force_contribution);
		}
		else
		{
			// Do nothing
			//node_iter->ClearAppliedForce();
		}
		//}
	}



	/*
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

		if( oxygenLevel > 0.99)
		{
			// Pressure is force per (cell diameter)^2
			// Force is to be applied towards the origin, radially
			c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(p_cell);
			c_vector<double, DIM> force_contribution;

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
		else
		{
			// Do nothing
			//node_iter->ClearAppliedForce();
		}



	}*/

}

template<unsigned DIM>
void ExternalPressureForceOnSurroundingSphere<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<mPressure>" << mPressure << "</mPressure> \n";

	// Call direct parent class
	AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ExternalPressureForceOnSurroundingSphere<1>;
template class ExternalPressureForceOnSurroundingSphere<2>;
template class ExternalPressureForceOnSurroundingSphere<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExternalPressureForceOnSurroundingSphere)

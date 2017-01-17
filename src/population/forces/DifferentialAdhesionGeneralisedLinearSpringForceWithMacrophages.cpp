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

#include "DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "NearMacrophageCellProperty.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages()
: GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>(),
  mHomotypicMacrophageSpringConstantMultiplier(1.0),
  mLabelledMacrophageSpringConstantMultiplier(1.0),
  mHomotypicLabelledSpringConstantMultiplier(1.0),
  mHomotypicUnlabelledSpringConstantMultiplier(1.0),
  mHeterotypicSpringConstantMultiplier(1.0)
  {
  }


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::VariableSpringConstantMultiplicationFactor(
		unsigned nodeAGlobalIndex,
		unsigned nodeBGlobalIndex,
		AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
		bool isCloserThanRestLength)
		{

	//if (isCloserThanRestLength)
	//{
	//	return 1.0;
	//}
	//else
	//{
		// Determine which (if any) of the cells corresponding to these nodes are labelled
		CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
		bool cell_A_is_labelled = p_cell_A->template HasCellProperty<NearMacrophageCellProperty>();

		CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);
		bool cell_B_is_labelled = p_cell_B->template HasCellProperty<NearMacrophageCellProperty>();

		if(cell_A_is_labelled or cell_B_is_labelled)
		{
			// What are the cell types? Note all macrophages should be labelled.
			bool cell_A_is_macrophage = p_cell_A->GetCellProliferativeType()->template IsType<MacrophageCellProliferativeType>();
			bool cell_B_is_macrophage = p_cell_B->GetCellProliferativeType()->template IsType<MacrophageCellProliferativeType>();
			if(cell_A_is_macrophage and cell_B_is_macrophage)
			{
				return mHomotypicMacrophageSpringConstantMultiplier;
			}
			if(cell_A_is_macrophage != cell_B_is_macrophage)
			{
				return mLabelledMacrophageSpringConstantMultiplier;
			}
			if(~cell_A_is_macrophage and ~cell_B_is_macrophage)
			{
				if(cell_A_is_labelled and cell_B_is_labelled)
				{
					return mHomotypicLabelledSpringConstantMultiplier;
				}
				else
				{
					return mHeterotypicSpringConstantMultiplier;
				}
			}
		}
		else
		{
			// If neither cell is labelled, we use HomotypicUnlabelledSpringConstantMultiplier
			return mHomotypicUnlabelledSpringConstantMultiplier;
		}
	//}

	// We should always return before reaching this line.
	NEVER_REACHED;
		}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::GetHomotypicMacrophageSpringConstantMultiplier()
{
	return mHomotypicMacrophageSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::SetHomotypicMacrophageSpringConstantMultiplier(double homotypicMacrophageSpringConstantMultiplier)
{
	assert(homotypicMacrophageSpringConstantMultiplier > 0.0);
	mHomotypicMacrophageSpringConstantMultiplier = homotypicMacrophageSpringConstantMultiplier;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::GetLabelledMacrophageSpringConstantMultiplier()
{
	return mLabelledMacrophageSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::SetLabelledMacrophageSpringConstantMultiplier(double labelledMacrophageSpringConstantMultiplier)
{
	assert(labelledMacrophageSpringConstantMultiplier > 0.0);
	mLabelledMacrophageSpringConstantMultiplier = labelledMacrophageSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::GetHomotypicLabelledSpringConstantMultiplier()
{
	return mHomotypicLabelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::SetHomotypicLabelledSpringConstantMultiplier(double homotypicLabelledSpringConstantMultiplier)
{
	assert(homotypicLabelledSpringConstantMultiplier > 0.0);
	mHomotypicLabelledSpringConstantMultiplier = homotypicLabelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::GetHomotypicUnlabelledSpringConstantMultiplier()
{
	return mHomotypicUnlabelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::SetHomotypicUnlabelledSpringConstantMultiplier(double homotypicUnlabelledSpringConstantMultiplier)
{
	assert(homotypicUnlabelledSpringConstantMultiplier > 0.0);
	mHomotypicUnlabelledSpringConstantMultiplier = homotypicUnlabelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::GetHeterotypicSpringConstantMultiplier()
{
	return mHeterotypicSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::SetHeterotypicSpringConstantMultiplier(double heterotypicSpringConstantMultiplier)
{
	assert(heterotypicSpringConstantMultiplier > 0.0);
	mHeterotypicSpringConstantMultiplier = heterotypicSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<HomotypicMacrophageSpringConstantMultiplier>" << mHomotypicMacrophageSpringConstantMultiplier << "</HomotypicMacrophageSpringConstantMultiplier>\n";
	*rParamsFile << "\t\t\t<LabelledMacrophageSpringConstantMultiplier>" << mLabelledMacrophageSpringConstantMultiplier << "</LabelledMacrophageSpringConstantMultiplier>\n";
	*rParamsFile << "\t\t\t<HomotypicLabelledSpringConstantMultiplier>" << mHomotypicLabelledSpringConstantMultiplier << "</HomotypicLabelledSpringConstantMultiplier>\n";
	*rParamsFile << "\t\t\t<HomotypicUnlabelledSpringConstantMultiplier>" << mHomotypicUnlabelledSpringConstantMultiplier << "</HomotypicUnlabelledSpringConstantMultiplier>\n";
	*rParamsFile << "\t\t\t<HeterotypicSpringConstantMultiplier>" << mHeterotypicSpringConstantMultiplier << "</HeterotypicSpringConstantMultiplier>\n";


	// Call direct parent class
	GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<1,1>;
template class DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<1,2>;
template class DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<2,2>;
template class DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<1,3>;
template class DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<2,3>;
template class DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages)

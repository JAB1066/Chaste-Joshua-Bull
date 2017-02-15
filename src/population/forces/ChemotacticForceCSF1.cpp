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

#include "ChemotacticForceCSF1.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

//#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"

template<unsigned DIM>
ChemotacticForceCSF1<DIM>::ChemotacticForceCSF1()
    : AbstractForce<DIM>(),
	sensitivity(1.0)
{
}

template<unsigned DIM>
ChemotacticForceCSF1<DIM>::~ChemotacticForceCSF1()
{
}

template<unsigned DIM>
double ChemotacticForceCSF1<DIM>::GetChemotacticForceMagnitude(const double concentration, const double concentrationGradientMagnitude)
{

    return sensitivity; // temporary force law - can be changed to something realistic
                          	  	  	  // without tests failing
}

template<unsigned DIM>
double ChemotacticForceCSF1<DIM>::GetChemotaxisSensitivity()
{
	return sensitivity;
}

template<unsigned DIM>
void ChemotacticForceCSF1<DIM>::SetChemotaxisSensitivity(double newSensitivity)
{
	sensitivity = newSensitivity;
}

void SetChemotaxisSensitivity(double newSensitivity);

template<unsigned DIM>
void ChemotacticForceCSF1<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
	// As CellwiseDataGradient is only defined for meshBasedCellPopulations, we have to calculate our own.
	// For each node, we get the corresponding element of the PDE mesh for CSF1. Then we look at the gradient on that element.

    //CellwiseDataGradient<DIM> gradients;
    //gradients.SetupGradients(rCellPopulation, "csf1");

	// Loop over cells, find which element of PDE mesh it is in, and determine gradient of PDE in that element
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	// Only macrophages move chemotactically
    	//boost::shared_ptr<AbstractCellProperty> p_prolifType = CellPropertyRegistry::Instance()->Get<MacrophageCellProliferativeType>();

        if (cell_iter->GetCellProliferativeType()->template IsType<MacrophageCellProliferativeType>())
        {
        	unsigned node_global_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        	// Find gradients in x, y, z directions at node
        	double xGrad = 0;
        	double yGrad = 0;
        	double zGrad = 0;



        	//gradient.clear(); // Ensure previous x, y, z gradients are cleared

        	xGrad = cell_iter->GetCellData()->GetItem("csf1_grad_x");
        	yGrad = cell_iter->GetCellData()->GetItem("csf1_grad_y");
        	if(DIM == 3)
        	{
            	zGrad = cell_iter->GetCellData()->GetItem("csf1_grad_z");
        	}

            //std::vector<double> gradient (3);
            c_vector<double,DIM> gradient;
        	gradient(0) = xGrad;
        	gradient(1) = yGrad;
        	if(DIM == 3)
        	{
        		gradient(2) = zGrad;
        	}


        	//gradient.push_back(yGrad);
        	//gradient.push_back(zGrad);

        	double magnitude_of_gradient = norm_2(gradient);
        	double nutrient_concentration = cell_iter->GetCellData()->GetItem("csf1");

        	// At the moment, magnitude of Chemotactic force is just the sensitivity
            double force_magnitude = GetChemotacticForceMagnitude(nutrient_concentration, magnitude_of_gradient);
            // force += chi * gradC/|gradC|
            if (magnitude_of_gradient > 0)
            {
                c_vector<double,DIM> force = force_magnitude*gradient;//(force_magnitude/magnitude_of_gradient)*gradient;
                rCellPopulation.GetNode(node_global_index)->AddAppliedForceContribution(force);
            }
            // else Fc=0
        }
    }
}

template<unsigned DIM>
void ChemotacticForceCSF1<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ChemotacticForceCSF1<1>;
template class ChemotacticForceCSF1<2>;
template class ChemotacticForceCSF1<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemotacticForceCSF1)

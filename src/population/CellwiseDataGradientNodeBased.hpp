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

#ifndef CELLWISEDATAGRADIENTNODEBASED_HPP_
#define CELLWISEDATAGRADIENTNODEBASED_HPP_

#include <NodeBasedCellPopulation.hpp>

/**
 *  A class for calculating the gradients of the CellwiseData, for Node Based cell populations.
 *
 *  Note - this is based on CellwiseDataGradient, but is designed to allow implementation of ChemotacticForceCSF1, which copies ChemotacticForce but for specific nutrient CSF1 and for node based cell populations.
 */
template<unsigned DIM>
class CellwiseDataGradientNodeBased
{
private:

    /**
     * The final gradients at the nodes
     */
    std::vector<c_vector<double, DIM> > mGradients;

public:

    /**
     * Compute the gradients at the nodes.
     *
     * We find the element of the PDE mesh containing this particular node, and then return the nutrient gradient there.
     *
     * @param rCellPopulation population on which to calculate gradients - must be instantiation of NodeBasedCellPopulation
     * @param rItemName is the name of the data from which to form the gradient (e.g. "oxygen"). (In this case, almost certainly CSF1)
     */
    void SetupGradients(AbstractCellPopulation<DIM>& rCellPopulation, const std::string& rItemName);

    /**
     * Get the gradient at a given node. Not set up for ghost nodes.
     *
     * @param nodeIndex
     *
     * @return the gradient at the node.
     */
    c_vector<double, DIM>& rGetGradient(unsigned nodeIndex);
};

#endif /*CELLWISEDATAGRADIENTNODEBASED_HPP_*/

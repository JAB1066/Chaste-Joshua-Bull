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

#ifndef TESTCUBOIDBOUNDARYCONDITIONS_HPP_
#define TESTCUBOIDBOUNDARYCONDITIONS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "DifferentiatedCellProliferativeType.hpp"
#include "CuboidPeriodicBoundaryCondition.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "DiffusionForceChooseD.hpp"
#include "NodeLocationWriter.hpp"
#include "OffLatticeSimulation.hpp"


#include "PetscSetupAndFinalize.hpp"

/**
 * This class contains tests for methods on classes inheriting from AbstractCellPopulationBoundaryCondition.
 */
class TestCuboidBoundaryConditions : public AbstractCellBasedTestSuite
{
public:

    void TestCuboidPeriodicBoundaryConditionWithNodeBasedCellPopulation() throw(Exception)
    {
    	// Generate Mesh:
    	// Make Vector
    	std::vector<Node<3>*> nodes;

    	// Add macrophage nodes in cube centred at origin, with edge size 10
    	for(unsigned nodeNum = 0; nodeNum < 1; nodeNum++)
    	{
    		double random_x = 10*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
    		double random_y = 10*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
    		double random_z = 10*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));

    		//nodes.push_back(new Node<3>(nodeNum, false, random_x, random_y, random_z));
    		nodes.push_back(new Node<3>(nodeNum, false, 0, 0, 0));

    	}

    	NodesOnlyMesh<3> mesh;
    	// Cut off length: 1.5 cell radii
    	mesh.ConstructNodesWithoutMesh(nodes, 1.5);

    	// Make cell pointers
    	std::vector<CellPtr> cells;

    	MAKE_PTR(WildTypeCellMutationState, p_state);
    	MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		// Now create macrophages
		for (unsigned i=0; i<nodes.size(); i++)
		{
			// Make Macrophage
			NoCellCycleModel* p_model = new NoCellCycleModel;
			p_model->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_macrophage_type);

			cells.push_back(p_cell);
		}

		// Make cell population
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<NodeLocationWriter>();

		OffLatticeSimulation<3> simulator(cell_population);

		// Create and set an output directory
		std::stringstream output_directory;
		output_directory << "AddingMacrophages/BoundaryConditions/PeriodicBCs/";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(300);

		ChastePoint<3> lower(-1,-1,-1);
		ChastePoint<3> upper(1,1,1);
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);

		simulator.Solve();

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
    }

};

#endif /*TESTCUBOIDBOUNDARYCONDITIONS_HPP_*/

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

#ifndef TESTARCHIVINGNODEBASEDSPHEROID_HPP_
#define TESTARCHIVINGNODEBASEDSPHEROID_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "DiffusionForce.hpp"
#include "DiffusionForceChooseD.hpp"


class TestArchivingNodeBasedSpheroid : public AbstractCellBasedTestSuite
{
public:

	void DontTestSaveArchiving() throw(Exception)
	{
		EXIT_IF_PARALLEL;

		std::vector<Node<3>*> nodes;

		nodes.push_back(new Node<3>(0,  false, 1, 1, 1));
		nodes.push_back(new Node<3>(1,  false, -1, -1, 1));
		nodes.push_back(new Node<3>(2,  false, 1, -1, -1));
		nodes.push_back(new Node<3>(3,  false, -1, 1, -1));

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);

		for (unsigned i=0; i<4; i++)
		{
			NoCellCycleModel* p_model = new NoCellCycleModel;
			CellPtr p_cell(new Cell(p_state, p_model));
			cells.push_back(p_cell);
		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		OffLatticeSimulation<3> simulator(cell_population);
		simulator.SetOutputDirectory("Temp/ArchivingTest");
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(0.5);

		// Add Brownian motion for all cells
		/* If both tests are run sequentially (save, then load), both tests pass for both of the below diffusion forces.
		 * However, if the archive is first created like this, and then you run only the load test, when using "DiffusionForce"
		 * it will load successfully, but when using "DiffusionForceChooseD" it will give a "what():  unregistered class" error
		 * */

		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		//MAKE_PTR(DiffusionForce<3>, p_diffusion_force);
		simulator.AddForce(p_diffusion_force);
		simulator.Solve();

		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);

		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void TestLoadArchiving() throw(Exception)
	{
		OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("Temp/ArchivingTest", 0.5);
		p_simulator->SetEndTime(1.5);
		p_simulator->Solve();
		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);
		delete p_simulator;
	}

};

#endif /*TESTARCHIVINGNODEBASEDSPHEROID_HPP_*/

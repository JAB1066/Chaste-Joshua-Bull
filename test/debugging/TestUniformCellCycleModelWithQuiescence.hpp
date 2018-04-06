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

#ifndef TESTUNIFORMCELLCYCLEMODELWITHQUIESCNCE_HPP_
#define TESTUNIFORMCELLCYCLEMODELWITHQUIESCNCE_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "UniformCellCycleModelWithQuiescence.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "NodesOnlyMesh.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "NodeLocationWriter.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "ApoptoticCellKiller.hpp"
//#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestUniformCellCycleModelWithQuiescence : public AbstractCellBasedTestSuite
{
public:


	void TestCellsProliferative() throw(Exception)
	{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes
		nodes.push_back(new Node<3>(0u,  false,  0.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(1u,  false,  2.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(2u,  false,  0.0, 2.0, 0.0));
		nodes.push_back(new Node<3>(3u,  false,  2.0, 2.0, 0.0));
		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

		// Create tumour cells manually
		for (unsigned i=0; i<4; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetMinCellCycleDuration(2.0);
			p_model->SetMaxCellCycleDuration(3.0);

			p_model->SetQuiescentConcentration(0.6);
			p_model->SetHypoxicConcentration(0.3);
			p_model->SetCriticalHypoxicDuration(1);

			CellPtr p_cell(new Cell(p_state, p_model));

			p_cell->SetCellProliferativeType(p_stem_type);

			p_cell->GetCellData()->SetItem("oxygen", 0.25*i);

			double birth_time = 0.0;
			p_cell->SetBirthTime(birth_time);

			cells.push_back(p_cell);
		}



		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<NodeLocationWriter>();

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.SetOutputDirectory("CellCycle/UniformModelWithQuiescence");//"TumourSpheroidSimulations/OffLattice3DTumourSpheroid");
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(10.0);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<3>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}



};

#endif /*TESTUNIFORMCELLCYCLEMODELWITHQUIESCNCE_HPP_*/

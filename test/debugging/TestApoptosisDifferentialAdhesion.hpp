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

#ifndef TESTAPOPTOSISDIFFERENTIALADHESION_HPP_
#define TESTAPOPTOSISDIFFERENTIALADHESION_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "DiffusionForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "UniformCellCycleModelWithQuiescence.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis.hpp"
#include "ApoptoticCellKiller.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "NodeLocationWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestArchiving : public AbstractCellBasedTestSuite
{
public:

	void TestApoptosisDifferentialAdhesion() throw(Exception)
	{
		const int DIM = 2;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		double width = 10;

		// Tumour Nodes
		unsigned nodeNum=0;
		for (double x=0; x<width+1; x++)
		{
			nodes.push_back(new Node<DIM>(nodeNum,  false,  x, 0));
			nodeNum++;
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

		// Create tumour cells manually
		for (unsigned i=0; i<nodeNum; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(DBL_MAX);
			p_model->SetMaxCellCycleDuration(DBL_MAX);
			p_model->SetHypoxicConcentration(0.5);
			p_model->SetQuiescentConcentration(0.5);
			p_model->SetCriticalHypoxicDuration(24);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);

			p_cell->GetCellData()->SetItem("oxygen", 1);
			if(i == 5)
			{
				p_cell->GetCellData()->SetItem("oxygen", 0);
			}
			p_cell->SetApoptosisTime(24); // Apoptosis time in hours - how long before a cell is removed?

			p_cell->SetBirthTime(0.0);
			cells.push_back(p_cell);

		}


		// Make cell population (1D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);

		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		std::stringstream output_directory;
		output_directory << "Debugging/DifferentialAdhesionApoptosis";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(60);
		simulator.SetEndTime(60);


		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(5.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
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

#endif /*TESTAPOPTOSISDIFFERENTIALADHESION_HPP_*/

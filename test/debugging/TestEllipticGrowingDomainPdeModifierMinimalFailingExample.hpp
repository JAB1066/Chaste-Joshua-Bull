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

#ifndef TESTELLIPTICGROWINGDOMAINPDEMODIFIERMINIMALFAILINGEXAMPLE_HPP_
#define TESTELLIPTICGROWINGDOMAINPDEMODIFIERMINIMALFAILINGEXAMPLE_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CheckpointArchiveTypes.hpp"
#include <boost/math/special_functions/bessel.hpp>

#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "UniformCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ReplicatableVector.hpp"
#include "OffLatticeSimulation.hpp"
#include "ApoptoticCellKiller.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"

#include "Debug.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * In this test suite we check the solution of the CellwisePdes against exact solutions.
 * In each case we are solving Laplacian U = f where f is constant in different regions.
 * We solve unit disc where the solutions are Bessel functions and logs.
 */
class TestEllipticGrowingDomainPdeModifier : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void dontTestNodeBasedSquareMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(20,20,0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 10.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),200u);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseEllipticPdeWithNodeOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA( p_cell_210->GetCellData()->GetItem("variable"), 0.4542, 1e-2); // Lower tolerance as slightly different meshes

        // Checking it doesn't change for this cell population
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.4476, 1e-4);

        // Clear memory
        delete p_mesh;
    }

    void TestEllipticGrowingDomainPdeModifierIn1d() throw(Exception)
        {
            // Create mesh
            std::vector<Node<1>*> nodes;
            nodes.push_back(new Node<1>(0, true,  0.0));
            nodes.push_back(new Node<1>(1, true, 1.0));
            nodes.push_back(new Node<1>(2, true,  2.0));
            nodes.push_back(new Node<1>(3, true,  3.0));
            nodes.push_back(new Node<1>(4, true,  4.0));
            nodes.push_back(new Node<1>(5, true,  5.0));
            nodes.push_back(new Node<1>(6, true,  6.0));
            nodes.push_back(new Node<1>(7, true,  7.0));
            nodes.push_back(new Node<1>(9, true,  8.0));

            NodesOnlyMesh<1> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);

            // Make cell pointers
            std::vector<CellPtr> cells;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(StemCellProliferativeType, p_stem_type);

            for (unsigned i=0; i<9; i++)
            		{
            			SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            			p_model->SetDimension(1);
            			p_model->SetHypoxicConcentration(0.5);
            			p_model->SetQuiescentConcentration(0.8);
            			p_model->SetCriticalHypoxicDuration(100);
            			p_model->SetStemCellG1Duration(0.1);
    					p_model->SetSDuration(0);
    					p_model->SetG2Duration(0);
    					p_model->SetMDuration(0);
            			CellPtr p_cell(new Cell(p_state, p_model));
            			p_cell->SetCellProliferativeType(p_stem_type);
            			p_cell->GetCellData()->SetItem("oxygen", 1.0);

            			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
            					(0.1);
            			p_cell->SetBirthTime(birthTime);
            			cells.push_back(p_cell);

            		}



            // Set up cell population
            NodeBasedCellPopulation<1> cell_population(mesh, cells);

            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(CellwiseSourceEllipticPde<1>, p_pde, (cell_population, -0.1));
            MAKE_PTR_ARGS(ConstBoundaryCondition<1>, p_bc, (1.0));;

            // Create a PDE modifier and set the name of the dependent variable in the PDE
            MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<1>, p_pde_modifier, (p_pde, p_bc, false));
            p_pde_modifier->SetDependentVariableName("oxygen");

            p_pde_modifier->SetOutputGradient(true);

            OffLatticeSimulation<1> simulator(cell_population);
            simulator.AddSimulationModifier(p_pde_modifier);
            std::stringstream output_directory;
            output_directory << "Debugging/1d_PDE";
            simulator.SetOutputDirectory(output_directory.str());
            simulator.SetSamplingTimestepMultiple(1);
            simulator.SetEndTime(10);

    		// Killer which removes apoptotic cells
    		MAKE_PTR_ARGS(ApoptoticCellKiller<1>, p_apoptosis_killer, (&cell_population));
    		simulator.AddCellKiller(p_apoptosis_killer);

    		MARK;
            simulator.Solve();
            MARK;
            // Clear memory
            for (unsigned i=0; i<nodes.size(); i++)
            {
                delete nodes[i];
            }
        }

    void dontTestEllipticGrowingDomainPdeModifierIn2d() throw(Exception)
    {
        // Create mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true,  0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(2, true,  2.0, 2.0));
        nodes.push_back(new Node<2>(3, true,  0.0, 4.0));
        nodes.push_back(new Node<2>(4, true,  4.0, 4.0));
        nodes.push_back(new Node<2>(5, true,  3.0, 3.0));
        nodes.push_back(new Node<2>(6, true,  3.0, 1.0));
        nodes.push_back(new Node<2>(7, true,  1.0, 3.0));
        nodes.push_back(new Node<2>(8, true,  1.0, 1.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned i=0; i<9; i++)
        		{
        			SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
        			p_model->SetDimension(2);
        			p_model->SetHypoxicConcentration(0.95);
        			p_model->SetQuiescentConcentration(0.95);
        			p_model->SetCriticalHypoxicDuration(0.2);
        			p_model->SetStemCellG1Duration(0.1);
					p_model->SetSDuration(0);
					p_model->SetG2Duration(0);
					p_model->SetMDuration(0);
        			CellPtr p_cell(new Cell(p_state, p_model));
        			p_cell->SetCellProliferativeType(p_stem_type);
        			p_cell->GetCellData()->SetItem("oxygen", 1.0);

        			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
        					(0.1);
        			p_cell->SetBirthTime(birthTime);
        			cells.push_back(p_cell);

        		}



        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));;

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        p_pde_modifier->SetOutputGradient(true);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.AddSimulationModifier(p_pde_modifier);
        std::stringstream output_directory;
        output_directory << "Debugging/pdeFailure";
        simulator.SetOutputDirectory(output_directory.str());
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(10);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

        simulator.Solve();
        // Clear memory
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void dontTestEllipticGrowingDomainPdeModifierIn3d() throw(Exception)
    {
        // Create a simple mesh
        TetrahedralMesh<3,3> temp_mesh;
        temp_mesh.ConstructCuboid(3,3,3);
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Set up cell population
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<3>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestEllipticGrowingDomainPdeModifierIn3d");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_62 = cell_population.GetCellUsingLocationIndex(13);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable"), 1.0, 1e-2);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable_grad_x"), 0.0, 1e-2);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable_grad_y"), 0.0, 1e-2);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable_grad_z"), 0.0, 1e-2);
    }
};

#endif /*TESTELLIPTICGROWINGDOMAINPDEMODIFIERMINIMALFAILINGEXAMPLE_HPP_*/

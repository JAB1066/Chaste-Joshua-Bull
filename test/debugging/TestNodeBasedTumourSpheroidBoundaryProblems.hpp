/*
	Joshua Bull, 18 October 2016, based on TestRunningNodeBasedSimulationsTutorial from CHASTE tutorials suite

 */

#ifndef TESTNODEBASEDTUMOURSPHEROIDBOUNDARYPROBLEMS_HPP_
#define TESTNODEBASEDTUMOURSPHEROIDBOUNDARYPROBLEMS_HPP_

#include <cxxtest/TestSuite.h>
//#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines a node-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"
// Write VTU files
#include "VoronoiDataWriter.hpp"
#include "RandomCellKiller.hpp"

// For test 3
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

#include "AveragedSourceEllipticPde.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticBoxDomainPdeModifierVariableTimestep.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"


class TestNodeBasedTumourSpheroidBoundaryProblems : public AbstractCellBasedTestSuite
{
public:

	void dontTestSpheroidWithOxygen() throw(Exception)
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<3>*> nodes;
        // Add some nodes

        nodes.push_back(new Node<3>(0u,  false,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4u,  false,  1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(5u,  false,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6u,  false,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(7u,  false,  1.0, 1.0, 1.0));
        NodesOnlyMesh<3> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        // Loop over nodes and create cells manually
        for (unsigned i=0; i<nodes.size(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(3);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1.0);

            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            p_model->SetStemCellG1Duration(4.0);
            p_model->SetHypoxicConcentration(0.1);
            p_model->SetQuiescentConcentration(0.3);
            p_model->SetCriticalHypoxicDuration(8);

            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  p_model->GetStemCellG1Duration()
                                  + p_model->GetSG2MDuration() );
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Make cell population (3D)
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files

        // Make PDE
        /*
         * With EllipticGrowingDomainPDE, BCs are correctly kept constant on cells on the outside. However, in the case the PDE remains constant throughout the domain.
         *
         */
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<3>, p_pde, (cell_population, -1.0)); // Cells should take in a lot of oxygen
        MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
        bool is_neumann_bc = false; // Dirichlet BCs

        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc));
        p_pde_modifier->SetDependentVariableName("oxygen");


        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.AddSimulationModifier(p_pde_modifier);
        simulator.SetOutputDirectory("NodeBasedTumourDebugging/GrowingDomain");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(30.0);

        /* Create force law */
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }


	void dontTestSpheroidWithOxygenInABox() throw(Exception)
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<3>*> nodes;
        // Add some nodes

        nodes.push_back(new Node<3>(0u,  false,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4u,  false,  1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(5u,  false,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6u,  false,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(7u,  false,  1.0, 1.0, 1.0));
        NodesOnlyMesh<3> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        // Loop over nodes and create cells manually
        for (unsigned i=0; i<nodes.size(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(3);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1.0);

            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            p_model->SetStemCellG1Duration(4.0);
            p_model->SetHypoxicConcentration(0.1);
            p_model->SetQuiescentConcentration(0.3);
            p_model->SetCriticalHypoxicDuration(8);

            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  p_model->GetStemCellG1Duration()
                                  + p_model->GetSG2MDuration() );
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Make cell population (3D)
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files


        /*
         * Simulation is identical to first test up to here
         * In this test, we create a box for PDE mesh and apply averagedSourceEllipticPde instead
         */

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        double cubeDomainDistanceToBoundary = 15;
        c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<3> lower(centroid(0)-cubeDomainDistanceToBoundary, centroid(1)-cubeDomainDistanceToBoundary, centroid(2)-cubeDomainDistanceToBoundary);
        ChastePoint<3> upper(centroid(0)+cubeDomainDistanceToBoundary, centroid(1)+cubeDomainDistanceToBoundary, centroid(2)+cubeDomainDistanceToBoundary);
        MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        // Make PDE
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -1.0)); //Again, cells take in lots of Oxygen
        MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
        bool is_neumann_bc = false; // Dirichlet BCs

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 5.0));
        p_pde_modifier->SetDependentVariableName("oxygen");
        p_pde_modifier->SetBcsOnBoxBoundary(true); // When this is set to false, BCs are applied to all cells. When set to true, cells on boundary decrease in oxygen concentration rather than stay at level of bcs.


        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
          * (this time with dimension 3) and set the output directory, output multiple and end time. */

         OffLatticeSimulation<3> simulator(cell_population);
         simulator.AddSimulationModifier(p_pde_modifier);
         simulator.SetOutputDirectory("NodeBasedTumourDebugging/InBox");
         simulator.SetSamplingTimestepMultiple(12);
         simulator.SetEndTime(30.0);

         /* Create force law */
         MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
         simulator.AddForce(p_force);

         simulator.Solve();

         /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
         for (unsigned i=0; i<nodes.size(); i++)
         {
             delete nodes[i];
         }
    }

	void TestVTKConcavehull() throw(Exception)
	    {
	        /** The next line is needed because we cannot currently run node based simulations in parallel. */
	        EXIT_IF_PARALLEL;

	        // Generate Mesh:
	        // Make Vector
	        std::vector<Node<2>*> nodes;
	        // Add some nodes

	        nodes.push_back(new Node<2>(0u,  false,  0.0, 0.0));
	        nodes.push_back(new Node<2>(1u,  false,  1.0, 0.0));
	        nodes.push_back(new Node<2>(2u,  false,  2.0, 0.0));
	        nodes.push_back(new Node<2>(3u,  false,  3.0, 0.0));
	        nodes.push_back(new Node<2>(4u,  false,  0.0, 1.0));
	        nodes.push_back(new Node<2>(5u,  false,  1.0, 1.0));
	        nodes.push_back(new Node<2>(6u,  false,  2.0, 1.0));
	        nodes.push_back(new Node<2>(7u,  false,  3.0, 1.0));
	        nodes.push_back(new Node<2>(8u,  false,  0.0, 0.0));
	        nodes.push_back(new Node<2>(9u,  false,  0.0, 1.0));
	        nodes.push_back(new Node<2>(10u,  false,  0.0, 2.0));
	        nodes.push_back(new Node<2>(11u,  false,  0.0, 3.0));

	        NodesOnlyMesh<2> mesh;
	        // Cut off length: 1.5 cell radii
	        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

	        // Make cell pointers
	        std::vector<CellPtr> cells;
	        MAKE_PTR(WildTypeCellMutationState, p_state);
	        MAKE_PTR(StemCellProliferativeType, p_stem_type);

	        // Loop over nodes and create cells manually
	        for (unsigned i=0; i<nodes.size(); i++)
	        {
	            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
	            p_model->SetDimension(2);

	            CellPtr p_cell(new Cell(p_state, p_model));
	            p_cell->SetCellProliferativeType(p_stem_type);
	            p_cell->GetCellData()->SetItem("oxygen", 1.0);

	            p_model->SetStemCellG1Duration(100.0);
	            p_model->SetTransitCellG1Duration(8.0);

	            p_model->SetStemCellG1Duration(4.0);
	            p_model->SetHypoxicConcentration(0.1);
	            p_model->SetQuiescentConcentration(0.3);
	            p_model->SetCriticalHypoxicDuration(8);

	            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
	                                 (  p_model->GetStemCellG1Duration()
	                                  + p_model->GetSG2MDuration() );
	            p_cell->SetBirthTime(birth_time);
	            cells.push_back(p_cell);
	        }

	        // Make cell population (3D)
	        NodeBasedCellPopulation<2> cell_population(mesh, cells);
	        cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files


	        /*
	         * Simulation is identical to first test up to here
	         * In this test, we create a box for PDE mesh and apply averagedSourceEllipticPde instead
	         */

	        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
	        double cubeDomainDistanceToBoundary = 15;
	        c_vector<double,2> centroid = cell_population.GetCentroidOfCellPopulation();
	        ChastePoint<2> lower(centroid(0)-cubeDomainDistanceToBoundary, centroid(1)-cubeDomainDistanceToBoundary, centroid(2)-cubeDomainDistanceToBoundary);
	        ChastePoint<2> upper(centroid(0)+cubeDomainDistanceToBoundary, centroid(1)+cubeDomainDistanceToBoundary, centroid(2)+cubeDomainDistanceToBoundary);
	        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

	        // Make PDE
	        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -1.0));
	        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));
	        bool is_neumann_bc = false; // Dirichlet BCs

	        // Create a PDE modifier and set the name of the dependent variable in the PDE
	        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<2>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 1.0));
	        p_pde_modifier->SetDependentVariableName("oxygen");
	        p_pde_modifier->SetBcsOnBoxBoundary(false); // When this is set to false, BCs are applied to all cells. When set to true, cells on boundary decrease in oxygen concentration rather than stay at level of bcs.


	        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
	          * (this time with dimension 3) and set the output directory, output multiple and end time. */

	         OffLatticeSimulation<2> simulator(cell_population);
	         simulator.AddSimulationModifier(p_pde_modifier);
	         simulator.SetOutputDirectory("NodeBasedTumourDebugging/VTK");
	         simulator.SetSamplingTimestepMultiple(1);
	         simulator.SetEndTime(2.0);

	         simulator.Solve();

	         /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
	         for (unsigned i=0; i<nodes.size(); i++)
	         {
	             delete nodes[i];
	         }
	    }

};


#endif /* TESTNODEBASEDTUMOURSPHEROIDBOUNDARYPROBLEMS_HPP_ */

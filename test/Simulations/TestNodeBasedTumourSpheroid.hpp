/*
	Joshua Bull, 18 October 2016, based on TestRunningNodeBasedSimulationsTutorial from CHASTE tutorials suite

 */

#ifndef TESTNODEBASEDTUMOURSPHEROID_HPP_
#define TESTNODEBASEDTUMOURSPHEROID_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
/* The remaining header files define classes that will be used in the cell population
 * simulation test. We encountered some of these header files in
 * UserTutorials/RunningMeshBasedSimulations. */
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"
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
/* Next, we define the test class.
 */
class TestNodeBasedTumourSpheroid : public AbstractCellBasedTestSuite
{
public:

	void TestSpheroid() throw(Exception)
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<3>*> nodes;
        // Add some nodes
        nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
        NodesOnlyMesh<3> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        // Make cell population (3D)
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TumourSpheroidSimulations/OffLattice3DTumourSpheroid");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(240.0);

        /* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        simulator.AddForce(p_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();


        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

};

#endif /* TESTNODEBASEDTUMOURSPHEROID_HPP_ */

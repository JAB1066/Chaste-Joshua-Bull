/*
	Joshua Bull, 18 October 2016, based on TestRunningNodeBasedSimulationsTutorial from CHASTE tutorials suite

 */

#ifndef TESTNODEBASEDTUMOURSPHEROID_HPP_
#define TESTNODEBASEDTUMOURSPHEROID_HPP_


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
#include "AveragedSourceParabolicPde.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"
#include "ApoptoticCellKiller.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic.hpp"
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"




/* Next, we define the test class.
 *
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
		CellsGenerator<UniformCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.SetOutputDirectory("Temp");//"TumourSpheroidSimulations/OffLattice3DTumourSpheroid");
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(60.0);

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


	void dontTestSpheroidWithRandomKiller() throw(Exception)
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
		CellsGenerator<UniformCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);


		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.SetOutputDirectory("TumourSpheroidSimulations/OffLattice3DTumourSpheroidWithCellDeath");
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(60.0);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		simulator.AddForce(p_force);

		// Now put in a random killer
		double probOfDeathPerHour;
		probOfDeathPerHour = 0.01;
		MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer, (&cell_population, probOfDeathPerHour));
		simulator.AddCellKiller(p_cell_killer);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
    						}


	void dontTestSpheroidWithOxygen() throw(Exception) // Temporarily broken
    						{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		unsigned nodeNum=0;
		double initialRadius = 3.0;
		for (double x=-initialRadius; x<initialRadius; x++)
		{
			for (double y=-initialRadius; y<initialRadius; y++)
			{
				for (double z=-initialRadius; z<initialRadius; z++)
				{
					if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
					{
						nodes.push_back(new Node<3>(nodeNum,  false,  x, y, z));
						nodeNum++;
					}
				}
			}
		}
		//        nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
		//        nodes.push_back(new Node<3>(1u,  false,  0.5, 1.0, 0.0));
		//        nodes.push_back(new Node<3>(2u,  false,  0.5, 0.0, 1.0));
		//        nodes.push_back(new Node<3>(3u,  false,  -0.5, 0.0, 0.5));
		//        nodes.push_back(new Node<3>(4u,  false,  0.0, 0.5, 0.0));
		//        nodes.push_back(new Node<3>(5u,  false,  0.0, -0.5, 0.5));
		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		//MAKE_PTR(WildTypeCellMutationState, p_state);
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);


		// Loop over nodes and create cells manually
		for (unsigned i=0; i<nodes.size(); i++)
		{
			SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
			p_model->SetDimension(3);
			//p_model->SetQuiescentConcentration(0.8); // Arbitrary
			//p_model->SetHypoxicConcentration(0.7); // Arbitrary

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 0.1);

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
		MAKE_PTR_ARGS(VolumeDependentAveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.03)); // Arbitrary - was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc));
		p_pde_modifier->SetDependentVariableName("oxygen");



		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("TumourSpheroidSimulations/OffLattice3DTumourSpheroidGrowingDomain");
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(240.0);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		simulator.AddForce(p_force);

		// Now put in a random killer
		double probOfDeathPerHour;
		probOfDeathPerHour = 0.01;
		MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer, (&cell_population, probOfDeathPerHour));
		simulator.AddCellKiller(p_cell_killer);

		/* To run the simulation, we call {{{Solve()}}}. */
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

		unsigned nodeNum=0;
		double initialRadius = 3.0;
		for (double x=-initialRadius; x<initialRadius; x++)
		{
			for (double y=-initialRadius; y<initialRadius; y++)
			{
				for (double z=-initialRadius; z<initialRadius; z++)
				{
					if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
					{
						nodes.push_back(new Node<3>(nodeNum,  false,  x, y, z));
						nodeNum++;
					}
				}
			}
		}
		//        nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
		//        nodes.push_back(new Node<3>(1u,  false,  0.5, 1.0, 0.0));
		//        nodes.push_back(new Node<3>(2u,  false,  0.5, 0.0, 1.0));
		//        nodes.push_back(new Node<3>(3u,  false,  -0.5, 0.0, 0.5));
		//        nodes.push_back(new Node<3>(4u,  false,  0.0, 0.5, 0.0));
		//        nodes.push_back(new Node<3>(5u,  false,  0.0, -0.5, 0.5));
		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		//MAKE_PTR(WildTypeCellMutationState, p_state);
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);


		// Loop over nodes and create cells manually
		for (unsigned i=0; i<nodes.size(); i++)
		{
			SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
			p_model->SetDimension(3);
			//p_model->SetQuiescentConcentration(0.8); // Arbitrary
			//p_model->SetHypoxicConcentration(0.7); // Arbitrary

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 0.1);

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
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		double cubeDomainDistanceToBoundary = 15;
		c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<3> lower(centroid(0)-cubeDomainDistanceToBoundary, centroid(1)-cubeDomainDistanceToBoundary, centroid(2)-cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(centroid(0)+cubeDomainDistanceToBoundary, centroid(1)+cubeDomainDistanceToBoundary, centroid(2)+cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 5.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("TumourSpheroidSimulations/OffLattice3DTumourSpheroidWithOxygenBoxBCs/PLOS13Parameters2");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(180.0);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Now put in a random killer
		double probOfDeathPerHour;
		probOfDeathPerHour = 0.01;
		MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer, (&cell_population, probOfDeathPerHour));
		simulator.AddCellKiller(p_cell_killer);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
    						}


	void dontTestSpheroidWithOxygenInABoxandCSF1() throw(Exception)
    			{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		unsigned nodeNum=0;
		double initialRadius = 7.0;
		for (double x=-initialRadius; x<initialRadius; x++)
		{
			for (double y=-initialRadius; y<initialRadius; y++)
			{
				for (double z=-initialRadius; z<initialRadius; z++)
				{
					if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
					{
						nodes.push_back(new Node<3>(nodeNum,  false,  x, y, z));
						nodeNum++;
					}
				}
			}
		}

		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		//MAKE_PTR(WildTypeCellMutationState, p_state);
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);


		// Loop over nodes and create cells manually
		for (unsigned i=0; i<nodes.size(); i++)
		{
			SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
			p_model->SetDimension(3);
			//p_model->SetQuiescentConcentration(0.8); // Arbitrary
			//p_model->SetHypoxicConcentration(0.7); // Arbitrary

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1.0);
			p_cell->GetCellData()->SetItem("csf1", 1.0);

			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);

			p_model->SetStemCellG1Duration(4.0);
			p_model->SetHypoxicConcentration(0.5);
			p_model->SetQuiescentConcentration(0.6);
			p_model->SetCriticalHypoxicDuration(8);//was 8

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(  p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files
		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
		cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();


		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Make PDE (CSF1)
		MAKE_PTR_ARGS(AveragedSourceParabolicPde<3>, p_pdeCSF, (cell_population, 1.0, 1.0, 0.03));//du/dt coefficient, diffusion coefficient, source coefficient
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bcCSF, (0.0));
		bool is_neumann_bcCSF = true; // Neumann BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		double cubeDomainDistanceToBoundary = 15;
		c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<3> lower(centroid(0)-cubeDomainDistanceToBoundary, centroid(1)-cubeDomainDistanceToBoundary, centroid(2)-cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(centroid(0)+cubeDomainDistanceToBoundary, centroid(1)+cubeDomainDistanceToBoundary, centroid(2)+cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<3>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
		p_pde_modifierCSF->SetDependentVariableName("csf1");
		p_pde_modifierCSF->SetBcsOnBoxBoundary(true); //was false


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		//simulator.AddSimulationModifier(p_pde_modifierCSF);
		simulator.SetOutputDirectory("TumourSpheroidSimulations/OffLattice3DTumourSpheroidWithOxygenBoxBCsandCSF/SteadyStateAttempt2");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(300.0);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Now put in a random killer
		double probOfDeathPerHour;
		probOfDeathPerHour = 0.01;
		MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer, (&cell_population, probOfDeathPerHour));
		simulator.AddCellKiller(p_cell_killer);

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


	void dontTestSpheroidWithOxygenInABoxandCSF1PauseCellsInAnyPhaseCellCycle() throw(Exception)
    			{
		// In previous tests, hypoxic cells continue their cell cycle until they reach the G1 phase. Here, we modify the cell cycle to cause them to freeze wherever they are in their cell cycle.
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		unsigned nodeNum=0;
		double initialRadius = 6.0;
		for (double x=-initialRadius; x<initialRadius; x++)
		{
			for (double y=-initialRadius; y<initialRadius; y++)
			{
				for (double z=-initialRadius; z<initialRadius; z++)
				{
					if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
					{
						nodes.push_back(new Node<3>(nodeNum,  false,  x, y, z));
						nodeNum++;
					}
				}
			}
		}

		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		//MAKE_PTR(WildTypeCellMutationState, p_state);
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);


		// Loop over nodes and create cells manually
		for (unsigned i=0; i<nodes.size(); i++)
		{
			SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic* p_model = new SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic;
			p_model->SetDimension(3);
			//p_model->SetQuiescentConcentration(0.8); // Arbitrary
			//p_model->SetHypoxicConcentration(0.7); // Arbitrary

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1.0);
			p_cell->GetCellData()->SetItem("csf1", 1.0);

			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);

			p_model->SetStemCellG1Duration(4.0);
			p_model->SetHypoxicConcentration(0.4);
			p_model->SetQuiescentConcentration(0.7);
			p_model->SetCriticalHypoxicDuration(6);//was 8

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(  p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files
		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
		cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();


		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Make PDE (CSF1)
		MAKE_PTR_ARGS(AveragedSourceParabolicPde<3>, p_pdeCSF, (cell_population, 1.0, 1.0, 0.03));//du/dt coefficient, diffusion coefficient, source coefficient
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bcCSF, (0.0));
		bool is_neumann_bcCSF = true; // Neumann BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		double cubeDomainDistanceToBoundary = 15;
		c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<3> lower(centroid(0)-cubeDomainDistanceToBoundary, centroid(1)-cubeDomainDistanceToBoundary, centroid(2)-cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(centroid(0)+cubeDomainDistanceToBoundary, centroid(1)+cubeDomainDistanceToBoundary, centroid(2)+cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<3>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
		p_pde_modifierCSF->SetDependentVariableName("csf1");
		p_pde_modifierCSF->SetBcsOnBoxBoundary(true); //was false


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		//simulator.AddSimulationModifier(p_pde_modifierCSF);
		simulator.SetOutputDirectory("TumourSpheroidSimulations/OffLattice3DTumourSpheroidWithOxygenBoxBCsandCSF/VaryCellCycle_2");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(180.0);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Now put in a random killer
		double probOfDeathPerHour;
		probOfDeathPerHour = 0.01;
		MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer, (&cell_population, probOfDeathPerHour));
		simulator.AddCellKiller(p_cell_killer);

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


	void dontTestSpheroidWithOxygenInABoxandCSF1Alarcon2004CellCycle() throw(Exception)
    			{
		// In previous tests, hypoxic cells continue their cell cycle until they reach the G1 phase. Here, we modify the cell cycle to cause them to freeze wherever they are in their cell cycle.
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		unsigned nodeNum=0;
		double initialRadius = 6.0;
		for (double x=-initialRadius; x<initialRadius; x++)
		{
			for (double y=-initialRadius; y<initialRadius; y++)
			{
				for (double z=-initialRadius; z<initialRadius; z++)
				{
					if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
					{
						nodes.push_back(new Node<3>(nodeNum,  false,  x, y, z));
						nodeNum++;
					}
				}
			}
		}

		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		//MAKE_PTR(WildTypeCellMutationState, p_state);
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);


		// Loop over nodes and create cells manually
		for (unsigned i=0; i<nodes.size(); i++)
		{
			Alarcon2004OxygenBasedCellCycleModel* p_model = new Alarcon2004OxygenBasedCellCycleModel;
			p_model->SetDimension(3);
			//p_model->SetQuiescentConcentration(0.8); // Arbitrary
			//p_model->SetHypoxicConcentration(0.7); // Arbitrary

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1.0);
			p_cell->GetCellData()->SetItem("csf1", 1.0);

			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);
			p_model->SetStemCellG1Duration(4.0);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(  p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files
		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
		cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();


		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Make PDE (CSF1)
		MAKE_PTR_ARGS(AveragedSourceParabolicPde<3>, p_pdeCSF, (cell_population, 1.0, 1.0, 0.03));//du/dt coefficient, diffusion coefficient, source coefficient
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bcCSF, (0.0));
		bool is_neumann_bcCSF = true; // Neumann BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		double cubeDomainDistanceToBoundary = 15;
		c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<3> lower(centroid(0)-cubeDomainDistanceToBoundary, centroid(1)-cubeDomainDistanceToBoundary, centroid(2)-cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(centroid(0)+cubeDomainDistanceToBoundary, centroid(1)+cubeDomainDistanceToBoundary, centroid(2)+cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<3>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
		p_pde_modifierCSF->SetDependentVariableName("csf1");
		p_pde_modifierCSF->SetBcsOnBoxBoundary(true); //was false


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		//simulator.AddSimulationModifier(p_pde_modifierCSF);
		simulator.SetOutputDirectory("TumourSpheroidSimulations/OffLattice3DTumourSpheroidWithOxygenBoxBCsandCSF/VaryCellCycleAlarcon");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(180.0);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Now put in a random killer
		double probOfDeathPerHour;
		probOfDeathPerHour = 0.01;
		MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer, (&cell_population, probOfDeathPerHour));
		simulator.AddCellKiller(p_cell_killer);

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


#endif /* TESTNODEBASEDTUMOURSPHEROID_HPP_ */

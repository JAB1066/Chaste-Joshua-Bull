#ifndef TESTFULLSPHEROIDSIMULATIONSWITHMACROPHAGES_HPP_
#define TESTFULLSPHEROIDSIMULATIONSWITHMACROPHAGES_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
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
#include "NodeLocationWriter.hpp"

#include "SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "ChemotacticForceCSF1.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

#include "AveragedSourceEllipticPde.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ApoptoticCellKiller.hpp"

#include "DiffusionForceChooseD.hpp"
#include "DiffusionForce.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "CuboidPeriodicBoundaryCondition.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



class TestFullSpheroidSimulationsWithMacrophages : public AbstractCellBasedTestSuite
{
public:

	void DontTestMacrophagesInfiltratingSpheroid() throw(Exception)
	{
		// In previous tests, hypoxic cells continue their cell cycle until they reach the G1 phase. Here, we modify the cell cycle to cause them to freeze wherever they are in their cell cycle.
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 25;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 9.0;
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

		unsigned firstMacrophageNode = nodeNum;

		// Macrophage Nodes - in a shell at fixed radius
		// Number of Macrophages remains constant throughout simulation
		unsigned numMacrophages=100;
		//double macrophageSphereRadius = 20.0;

		std::vector<double> xVec;
		std::vector<double> yVec;
		std::vector<double> zVec;
		unsigned macCount=0;
		double random_x;
		double random_y;
		double zCoord;
		/*
		while(i < numMacrophages)
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			random_x = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			random_y = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			xVec.push_back(random_x);
			yVec.push_back(random_y);

			// Choose sign of z coord randomly
			double randNum = RandomNumberGenerator::Instance()->ranf();
			if(randNum < 0.5)
			{
				zCoord = sqrt(abs(pow(macrophageSphereRadius,2) - pow(random_x,2) - pow(random_y,2)));
			}
			else
			{
				zCoord = -sqrt(abs(pow(macrophageSphereRadius,2) - pow(random_x,2) - pow(random_y,2)));
			}
			zVec.push_back(zCoord);
			i++;
		}
		 */
		// Add macrophages instead at random points in cube outside of spheroid radius
		while(macCount<numMacrophages)
		{
			random_x = cubeDomainDistanceToBoundary*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			random_y = cubeDomainDistanceToBoundary*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			zCoord = cubeDomainDistanceToBoundary*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			if((pow(random_x,2)+pow(random_y,2)+pow(zCoord,2)) > pow(initialRadius,2))
			{
				xVec.push_back(random_x);
				yVec.push_back(random_y);
				zVec.push_back(zCoord);

				macCount++;
			}
		}

		for(unsigned i=0; i<xVec.size();i++)
		{
			nodes.push_back(new Node<3>(nodeNum,  false,  xVec[i], yVec[i], zVec[i]));
			nodeNum++;
		}

		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		//MAKE_PTR(WildTypeCellMutationState, p_state);
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);


		// Create tumour cells manually
		for (unsigned i=0; i<firstMacrophageNode; i++)
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
			p_model->SetHypoxicConcentration(0.3);
			p_model->SetQuiescentConcentration(0.4);
			p_model->SetCriticalHypoxicDuration(6);//was 8

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(  p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
		// Now create macrophages
		for (unsigned i=firstMacrophageNode; i<nodes.size(); i++)
		{
			// Make Macrophage
			NoCellCycleModel* p_model2 = new NoCellCycleModel;
			p_model2->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model2));
			p_cell->SetCellProliferativeType(p_macrophage_type);
			p_cell->GetCellData()->SetItem("csf1", 1.0);

			cells.push_back(p_cell);

		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		// Write summary statistic files
		//cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
		//cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
		//cell_population.AddCellWriter<CellMutationStatesWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();

		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Make PDE (CSF1)
		//MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pdeCSF, (cell_population, 0.03));
		MAKE_PTR_ARGS(AveragedSourceParabolicPde<3>, p_pdeCSF, (cell_population, 1.0, 1.0, 0.03));//du/dt coefficient, diffusion coefficient, source coefficient
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bcCSF, (0.0));
		bool is_neumann_bcCSF = true; // Neumann BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<3> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		//MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<3>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 2.0));
		p_pde_modifierCSF->SetDependentVariableName("csf1");
		p_pde_modifierCSF->SetOutputGradient(true); // Must be set to true for CSF chemotaxis to work
		p_pde_modifierCSF->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.AddSimulationModifier(p_pde_modifierCSF);
		simulator.SetOutputDirectory("FullSimulations/InitialTesting/6");
		simulator.SetSamplingTimestepMultiple(36);//12);
		simulator.SetEndTime(100.0);

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.2);
		simulator.AddForce(p_diffusion_force);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<3>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		// Add Chemotaxis for macrophages
		MAKE_PTR(ChemotacticForceCSF1<3>, p_chemotactic_force);
		p_chemotactic_force->SetChemotaxisSensitivity(0.2);
		simulator.AddForce(p_chemotactic_force);


		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void TestSpheroidGrowthCreateBenchmarkTumour() throw(Exception)
    															{
		/* This is a simulation designed to provide groundwork to reproducing Dorie et al.
		 * "Migration and Internalization of Cells and Polystyrene Microspheres in Tumour Cell Spheroids"
		 * Experimental Cell Research 141 (1982) 201-209
		 *
		 * In this case, inert polystyrene beads infiltrate a tumour spheroid where the only forces at work are proliferation and brownian motion
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 25;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 2;//13.5;
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

		unsigned firstMacrophageNode = nodeNum;

		// Macrophage Nodes - in a shell at fixed radius
		// Number of Macrophages remains constant throughout simulation
		unsigned numMacrophages=0;

		std::vector<double> xVec;
		std::vector<double> yVec;
		std::vector<double> zVec;
		unsigned macCount=0;
		double random_x;
		double random_y;
		double zCoord;
		/*
		while(i < numMacrophages)
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			random_x = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			random_y = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			xVec.push_back(random_x);
			yVec.push_back(random_y);

			// Choose sign of z coord randomly
			double randNum = RandomNumberGenerator::Instance()->ranf();
			if(randNum < 0.5)
			{
				zCoord = sqrt(abs(pow(macrophageSphereRadius,2) - pow(random_x,2) - pow(random_y,2)));
			}
			else
			{
				zCoord = -sqrt(abs(pow(macrophageSphereRadius,2) - pow(random_x,2) - pow(random_y,2)));
			}
			zVec.push_back(zCoord);
			i++;
		}
		 */
		// Add macrophages instead at random points in cube outside of spheroid radius
		while(macCount<numMacrophages)
		{
			random_x = cubeDomainDistanceToBoundary*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			random_y = cubeDomainDistanceToBoundary*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			zCoord = cubeDomainDistanceToBoundary*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			if((pow(random_x,2)+pow(random_y,2)+pow(zCoord,2)) > pow(initialRadius,2))
			{
				xVec.push_back(random_x);
				yVec.push_back(random_y);
				zVec.push_back(zCoord);
				macCount++;
			}
		}

		for(unsigned i=0; i<xVec.size();i++)
		{
			nodes.push_back(new Node<3>(nodeNum,  false,  xVec[i], yVec[i], zVec[i]));
			nodeNum++;
		}

		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		//MAKE_PTR(WildTypeCellMutationState, p_state);
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);


		// Create tumour cells manually
		for (unsigned i=0; i<firstMacrophageNode; i++)
		{
//			SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic* p_model = new SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic;
			//TODO Archiving problem
			NoCellCycleModel* p_model = new NoCellCycleModel;
//			p_model->SetDimension(3);
//			//p_model->SetQuiescentConcentration(0.8); // Arbitrary
//			//p_model->SetHypoxicConcentration(0.7); // Arbitrary

			CellPtr p_cell(new Cell(p_state, p_model));
//			p_cell->SetCellProliferativeType(p_stem_type);
//			p_cell->GetCellData()->SetItem("oxygen", 1.0);
//
//			p_model->SetStemCellG1Duration(8.0);
//			p_model->SetTransitCellG1Duration(8.0);
//
//			p_model->SetStemCellG1Duration(4.0);
//			p_model->SetHypoxicConcentration(0.2);
//			p_model->SetQuiescentConcentration(0.4);
//			p_model->SetCriticalHypoxicDuration(6);//was 8
//
//			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
//					(  p_model->GetStemCellG1Duration()
//							+ p_model->GetSG2MDuration() );
//			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
		// Now create macrophages
		for (unsigned i=firstMacrophageNode; i<nodes.size(); i++)
		{
			// Make Macrophage
			NoCellCycleModel* p_model2 = new NoCellCycleModel;
			p_model2->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model2));
			p_cell->SetCellProliferativeType(p_macrophage_type);
			p_cell->GetCellData()->SetItem("csf1", 1.0);

			cells.push_back(p_cell);

		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();

		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<3> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		//MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 1.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
		//PRINT_VARIABLE(p_pde_modifier->GetFeMesh()->GetNumNodes());


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		//simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("FullSimulations/Dorie1982/SpheroidOnly/one_hour_archive_test");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(0.1);

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		//simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);
		//TODO Archiving Problem

		// Add Brownian motion for all cells
		//MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		MAKE_PTR(DiffusionForce<3>, p_diffusion_force);
		//p_diffusion_force->SetDiffusionScalingConstant(0.2);
		simulator.AddForce(p_diffusion_force);
		//TODO Archiving Problem in "choose D"

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<3>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();
		//PRINT_VARIABLE(p_pde_modifier->GetFeMesh()->GetNumNodes());

		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);
		//std::cout << "Simulation Modifiers size = " << simulator.GetSimulationModifiers()->size() << std::endl;

//		MARK;
//		OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("FullSimulations/Dorie1982/SpheroidOnly/one_hour_archive_test", 0.2);
//		MARK;
//		p_simulator->SetEndTime(0.4);
//		p_simulator->Solve();
//		MARK;
//		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);
//		MARK;
//		delete p_simulator;

//		MARK;
//		OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("FullSimulations/Dorie1982/SpheroidOnly/one_hour_archive_test", 0.2);
//		std::cout << "Simulation Modifiers size = " << p_simulator->GetSimulationModifiers()->size() << std::endl;
//		MARK;
//		p_simulator->SetEndTime(0.4);
//		MARK;
//		PRINT_VARIABLE(p_pde_modifier->GetFeMesh()->GetNumNodes());
//		//EllipticBoxDomainPdeModifier<3> tempModifier = p_simulator->GetSimulationModifiers()[0];
//		//PRINT_VARIABLE(tempModifier->GetFeMesh()->GetNumNodes());
//		MARK;
//		p_simulator->Solve();
//		MARK;
//		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);
//		MARK;
//		delete p_simulator;

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
    															}

	void TestSpheroidGrowthLoadBenchmarkTumour() throw(Exception)
	    															{
		/* This is a simulation designed to provide groundwork to reproducing Dorie et al.
		 * "Migration and Internalization of Cells and Polystyrene Microspheres in Tumour Cell Spheroids"
		 * Experimental Cell Research 141 (1982) 201-209
		 *
		 * In this case, inert polystyrene beads infiltrate a tumour spheroid where the only forces at work are proliferation and brownian motion
		 */


		MARK;
		OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("FullSimulations/Dorie1982/SpheroidOnly/one_hour_archive_test", 0.1);
		std::cout << "Simulation Modifiers size = " << p_simulator->GetSimulationModifiers()->size() << std::endl;
		MARK;

		//p_simulator->GetSimulationModifiers()->GenerateFeMesh;
//		double cubeDomainDistanceToBoundary = 25;
//		ChastePoint<3> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
//		ChastePoint<3> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
//		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));
//
//		p_simulator->GetSimulationModifiers()->GenerateFeMesh(p_cuboid,1.0);

		//p_simulator->GetSimulationModifiers()[0];


		p_simulator->SetEndTime(0.4);
		p_simulator->Solve();
		MARK;
		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);
		MARK;
		delete p_simulator;



	    															}

};


#endif /* TESTFULLSPHEROIDSIMULATIONSWITHMACROPHAGES_HPP_ */

#ifndef TESTSPHEROIDADDINGMACROPHAGES_HPP_
#define TESTSPHEROIDADDINGMACROPHAGES_HPP_


#include <cxxtest/TestSuite.h>
//#include "CheckpointArchiveTypes.hpp"
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
#include "VoronoiDataWriter.hpp"
#include "RandomCellKiller.hpp"

// For test 3
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "ChemotacticForceCSF1.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

#include "AveragedSourceEllipticPde.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"
#include "ApoptoticCellKiller.hpp"

#include "DiffusionForce.hpp"
#include "DiffusionForceChooseD.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages.hpp"
#include "MacrophageProximityLabellingModifier.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>


/**
 * Functional Boundary Condition, setting BC = bcConc if z = distance
 */
double bc_func(const ChastePoint<3>& p)
{
	double bcConc = 1.0;
	double distance = 10;
	double value = 0.0;

	if(p[2] == distance)
	{
		value = bcConc;
	}

	return value;
}

double bc_func2(const ChastePoint<2>& p)
{
	double bcConc = 1.0;
	double distance = 15;
	double value = 0.0;

	if(p[1] == distance)
	{
		value = bcConc;
	}

	return value;
}

class TestSpheroidAddingMacrophages : public AbstractCellBasedTestSuite
{
public:

	void dontTestAddSecondCellTypeToSpheroid() throw(Exception)
	{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		unsigned nodeNum=0;
		double initialRadius = 5.0;
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

		//boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);


		// Loop over nodes and create cells manually
		for (unsigned i=0; i<nodes.size(); i++)
		{
			SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
			p_model->SetDimension(3);


			CellPtr p_cell(new Cell(p_state, p_model));
			if(i % 3 == 0)
			{
				p_cell->SetCellProliferativeType(p_differentiated_type);
			}
			else
			{
				p_cell->SetCellProliferativeType(p_stem_type);
			}

			p_cell->GetCellData()->SetItem("oxygen", 1.0);
			p_cell->GetCellData()->SetItem("csf1", 1.0);

			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);

			p_model->SetStemCellG1Duration(4.0);
			p_model->SetHypoxicConcentration(0.5);
			p_model->SetQuiescentConcentration(0.7);
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
		simulator.SetOutputDirectory("AddingMacrophages/TwoCellTypes");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(90.0);

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

	void dontTestMacrophageRespondingToCSF1() throw(Exception)
															{
		// This test simulates a macrophage moving up a chemotactic gradient in 3D.
		EXIT_IF_PARALLEL;

		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add macrophage node

		nodes.push_back(new Node<3>(0,  false,  0, 0, 0));
		nodes.push_back(new Node<3>(1,  false,  -5, 0, 0));
		nodes.push_back(new Node<3>(2,  false,  5, 0, 0));


		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;

		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		//boost::shared_ptr<AbstractCellProperty> p_prolifType = CellPropertyRegistry::Instance()->Get<MacrophageCellProliferativeType>();

		// Loop over nodes and create cells manually
		for (unsigned i=0; i<nodes.size(); i++)
		{
			NoCellCycleModel* p_model = new NoCellCycleModel;
			p_model->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_macrophage_type);

			p_cell->GetCellData()->SetItem("csf1", 1.0);
			cells.push_back(p_cell);
		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		double cubeDomainDistanceToBoundary = 10;
		//c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
		c_vector<double,3> centroid;
		centroid(0)=0;
		centroid(1)=0;
		centroid(2)=0; // IF YOU CHANGE DISTANCE, ALSO CHANGE BC FUNCTION!!!!!
		ChastePoint<3> lower(centroid(0)-cubeDomainDistanceToBoundary, centroid(1)-cubeDomainDistanceToBoundary, centroid(2)-cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(centroid(0)+cubeDomainDistanceToBoundary, centroid(1)+cubeDomainDistanceToBoundary, centroid(2)+cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));


		// Make PDE
		//UniformSourceEllipticPde<3> pde(0);
		MAKE_PTR_ARGS(UniformSourceEllipticPde<3>, p_pde, (0.0));
		MAKE_PTR_ARGS(FunctionalBoundaryCondition<3>, p_functional_bc, (&bc_func));
		//MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_functional_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("csf1");
		p_pde_modifier->SetOutputGradient(true);
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		// Make Chemotactic Force for CSF1
		MAKE_PTR(ChemotacticForceCSF1<3>, p_chemotactic_force);

		//MAKE_PTR(DiffusionForce<3>, p_diffusion_force);

		//MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		//p_force->SetMeinekeSpringStiffness(30.0);
		//p_force->SetCutOffLength(1.5);
		//simulator.AddForce(p_force);

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		//simulator.AddForce(p_diffusion_force);
		simulator.AddForce(p_chemotactic_force);
		simulator.SetOutputDirectory("AddingMacrophages/WithChemotaxis");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(300.0);


		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
															}

	void dontTestMacrophagePassingThroughMonolayer() throw(Exception)
														{
		// This test simulates a monolayer of (non-proliferating) tumour cells with a macrophage at one end and a source of CSF at the other.
		//The macrophage should move up the chemotactic gradient.

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<2>*> nodes;

		// Add some nodes
		//Construct monolayer
		unsigned nodeNum=0;
		double halfWidth = 10.0;
		double halfDepth = 2.5;
		for (double x=-halfWidth; x<halfWidth; x++)
		{
			for (double y=-halfDepth; y<halfDepth; y++)
			{
				nodes.push_back(new Node<2>(nodeNum,  false,  x, y));
				nodeNum++;
			}
		}

		// Add macrophage node
		nodes.push_back(new Node<2>(nodeNum, false, 0, -halfDepth*2));

		NodesOnlyMesh<2> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;

		//boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);


		// Loop over nodes and create cells manually
		for (unsigned i=0; i<(nodes.size()-1); i++)
		{
			NoCellCycleModel* p_model = new NoCellCycleModel;
			p_model->SetDimension(2);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_differentiated_type);
			p_cell->GetCellData()->SetItem("csf1", 1.0);

			cells.push_back(p_cell);
		}

		// Make Macrophage
		NoCellCycleModel* p_model = new NoCellCycleModel;
		p_model->SetDimension(2);

		CellPtr p_cell(new Cell(p_state, p_model));
		p_cell->SetCellProliferativeType(p_macrophage_type);
		p_cell->GetCellData()->SetItem("csf1", 1.0);

		cells.push_back(p_cell);



		// Make cell population (2D)
		NodeBasedCellPopulation<2> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files
		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		double cubeDomainDistanceToBoundary = 15;
		ChastePoint<2> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<2> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

		// Make PDE
		//UniformSourceEllipticPde<3> pde(0);
		MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (0.0));
		MAKE_PTR_ARGS(FunctionalBoundaryCondition<2>, p_functional_bc, (&bc_func2));
		//MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_functional_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("csf1");
		p_pde_modifier->SetOutputGradient(true);
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 2) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<2> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("AddingMacrophages/ChemotaxisThroughMonolayer");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(150.0);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		/*
		 *
		 MAKE_PTR(DiffusionForce<2>, p_diffusion_force);
		//p_diffusion_force->SetViscosity();
		p_diffusion_force->SetAbsoluteTemperature(1);
		simulator.AddForce(p_diffusion_force);
		 *
		 */

		// Chemotaxis
		MAKE_PTR(ChemotacticForceCSF1<2>, p_chemotactic_force);
		p_chemotactic_force->SetChemotaxisSensitivity(1.0);
		simulator.AddForce(p_chemotactic_force);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
														}

	void dontTestMacrophagePassingThroughMonolayerDifferentialAdhesion() throw(Exception)
														{
		// This test simulates a monolayer of (non-proliferating) tumour cells with a macrophage at one end and a source of CSF at the other.
		//The macrophage should move up the chemotactic gradient. We include differential adhesion indicating weakened intercellular adhesion in the presence of macrophages

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<2>*> nodes;

		// Add some nodes
		//Construct monolayer
		unsigned nodeNum=0;
		double halfWidth = 10.0;
		double halfDepth = 2.5;
		for (double x=-halfWidth; x<halfWidth; x++)
		{
			for (double y=-halfDepth; y<halfDepth; y++)
			{
				nodes.push_back(new Node<2>(nodeNum,  false,  x, y));
				nodeNum++;
			}
		}

		// Add macrophage node
		nodes.push_back(new Node<2>(nodeNum, false, 0, -halfDepth*2));

		NodesOnlyMesh<2> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;

		//boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);


		// Loop over nodes and create cells manually
		for (unsigned i=0; i<(nodes.size()-1); i++)
		{
			NoCellCycleModel* p_model = new NoCellCycleModel;
			p_model->SetDimension(2);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_differentiated_type);
			p_cell->GetCellData()->SetItem("csf1", 1.0);

			cells.push_back(p_cell);
		}

		// Make Macrophage
		NoCellCycleModel* p_model = new NoCellCycleModel;
		p_model->SetDimension(2);

		CellPtr p_cell(new Cell(p_state, p_model));
		p_cell->SetCellProliferativeType(p_macrophage_type);
		p_cell->GetCellData()->SetItem("csf1", 1.0);

		cells.push_back(p_cell);


		// Make cell population (2D)
		NodeBasedCellPopulation<2> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files
		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		double cubeDomainDistanceToBoundary = 15;
		ChastePoint<2> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<2> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

		// Make PDE
		//UniformSourceEllipticPde<3> pde(0);
		MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (0.0));
		MAKE_PTR_ARGS(FunctionalBoundaryCondition<2>, p_functional_bc, (&bc_func2));
		//MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_functional_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("csf1");
		p_pde_modifier->SetOutputGradient(true);
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 2) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<2> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("AddingMacrophages/ChemotaxisThroughMonolayerDifferentialAdhesion");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(150.0);


		/* Again we create a force law (this time with dimension 2), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<2>, p_force);
		// Anything within the set distance of a macrophage has half the spring force of a usual tumour-tumour spring.
		//p_force->SetHomotypicMacrophageSpringConstantMultiplier(1.0);
		p_force->SetLabelledMacrophageSpringConstantMultiplier(1.0);
		p_force->SetHomotypicLabelledSpringConstantMultiplier(0.5);
		p_force->SetHeterotypicSpringConstantMultiplier(1);

		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);

		// Before adding force law, we need the helper modifier MacrophageProximityLabellingModifier
		MAKE_PTR(MacrophageProximityLabellingModifier<2>, p_macProx);
		p_macProx->SetMacrophageProximityLabellingRadius(5.0);

		simulator.AddSimulationModifier(p_macProx);
		simulator.AddForce(p_force);


		MAKE_PTR(DiffusionForce<2>, p_diffusion_force);
		//p_diffusion_force->SetViscosity();
		p_diffusion_force->SetAbsoluteTemperature(1);
		simulator.AddForce(p_diffusion_force);


		// Chemotaxis
		MAKE_PTR(ChemotacticForceCSF1<2>, p_chemotactic_force);
		p_chemotactic_force->SetChemotaxisSensitivity(2.0);
		simulator.AddForce(p_chemotactic_force);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
														}

	void dontTestMacrophagesInfiltratingSpheroid() throw(Exception)
    													{
		// In previous tests, hypoxic cells continue their cell cycle until they reach the G1 phase. Here, we modify the cell cycle to cause them to freeze wherever they are in their cell cycle.
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 8.0;
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
		unsigned numMacrophages=100;
		double macrophageSphereRadius = 10.0;

		std::vector<double> xVec;
		std::vector<double> yVec;
		std::vector<double> zVec;
		unsigned i=0;
		double random_x;
		double random_y;
		double zCoord;
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
			p_model->SetHypoxicConcentration(0.4);
			p_model->SetQuiescentConcentration(0.7);
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
		cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files
		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
		cell_population.AddCellWriter<CellMutationStatesWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();

		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Make PDE (CSF1)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pdeCSF, (cell_population, 0.03));
		//MAKE_PTR_ARGS(AveragedSourceParabolicPde<3>, p_pdeCSF, (cell_population, 1.0, 1.0, 0.03));//du/dt coefficient, diffusion coefficient, source coefficient
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bcCSF, (0.0));
		bool is_neumann_bcCSF = true; // Neumann BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		double cubeDomainDistanceToBoundary = 20;
		c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<3> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		//MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<3>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 2.0));
		p_pde_modifierCSF->SetDependentVariableName("csf1");
		p_pde_modifierCSF->SetOutputGradient(true);
		p_pde_modifierCSF->SetBcsOnBoxBoundary(true); //was false




		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.AddSimulationModifier(p_pde_modifierCSF);
		simulator.SetOutputDirectory("AddingMacrophages/InfiltrationIntoSpheroid");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(30.0);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		MAKE_PTR(DiffusionForce<3>, p_diffusion_force);
		//p_diffusion_force->SetViscosity();
		p_diffusion_force->SetAbsoluteTemperature(10);
		simulator.AddForce(p_diffusion_force);

		// Add Chemotaxis
		MAKE_PTR(ChemotacticForceCSF1<3>, p_chemotactic_force);
		simulator.AddForce(p_chemotactic_force);

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

	void TestMacrophageBrownianMotion() throw(Exception)
    													{
		// Test macrophage random motion and log location after time t
		/* We vary parameters:
		 *
		 * diffusionCoefficient: coefficient for random macrophage motion law (Fickian)
		 * simulationDuration: length of simulation in hours
		 *
		 */
		const int dimensions = 2;
		double diffusionCoefficient = 1.0; // D
		double simulationDuration = 120;
		std::string id_string = "1";

		// Generate Mesh:
		// Make Vector
		std::vector<Node<dimensions>*> nodes;

		// Add macrophage node at origin
		unsigned nodeNum=0;

		// Add macrophage node
		switch(dimensions)
		{
		case 1:
			nodes.push_back(new Node<dimensions>(nodeNum, false, 0));
			break;

		case 2:
			nodes.push_back(new Node<dimensions>(nodeNum, false, 0, 0));
			break;

		case 3:
			nodes.push_back(new Node<dimensions>(nodeNum, false, 0, 0, 0));
			break;
		}


		NodesOnlyMesh<dimensions> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;

		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		// Make Macrophage
		NoCellCycleModel* p_model = new NoCellCycleModel;
		p_model->SetDimension(dimensions);

		CellPtr p_cell(new Cell(p_state, p_model));
		p_cell->SetCellProliferativeType(p_macrophage_type);

		cells.push_back(p_cell);

		// Make cell population (2D)
		NodeBasedCellPopulation<dimensions> cell_population(mesh, cells);
		//cell_population.AddPopulationWriter<NodeLocationWriter>();

		cell_population.SetDampingConstantNormal(1.0);
		cell_population.SetDampingConstantMutant(1.0);


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 2) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<dimensions> simulator(cell_population);

		// Create and set an output directory that is different for each simulation
		std::stringstream output_directory;
		output_directory << "AddingMacrophages/DiffusionImplementationTesting/SingleMacrophage_DistanceInSetTime_OriginalRule/" << id_string;
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(simulationDuration);

		MAKE_PTR(DiffusionForceChooseD<dimensions>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(diffusionCoefficient);
		simulator.AddForce(p_diffusion_force);


		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();

		OutputFileHandler results_handler(output_directory.str(), false);
		out_stream results_file = results_handler.OpenOutputFile("results.completiontimeandlocation");

		Node<dimensions>* p_MacNode = cell_population.GetNode(nodeNum);
		const c_vector<double, dimensions>& r_MacNode_location = p_MacNode->rGetLocation();



		// Output summary statistics to results file
		(*results_file) << "ID" << ","
						<< "Simulation error" << ",";
		switch(dimensions)
		{
		case(1):
				(*results_file) << "Macrophage x coordinate at end" << ",";
		break;
		case(2):
				(*results_file) << "Macrophage x coordinate at end" << ","
								<< "Macrophage y coordinate at end" << ",";
		break;
		case(3):
				(*results_file) << "Macrophage x coordinate at end" << ","
								<< "Macrophage y coordinate at end" << ","
								<< "Macrophage z coordinate at end" << ",";
		break;
		}

		(*results_file) << "Total Macrophage displacement" << ","
						<< "Simulation time at end" << std::endl;

		double displacement;

		(*results_file) << id_string << ","
				<< "0" << ",";
		switch(dimensions)
		{
		case(1):
				(*results_file) << boost::lexical_cast<std::string>(r_MacNode_location[0]) << ",";
				displacement = sqrt(r_MacNode_location[0]*r_MacNode_location[0]);
		break;
		case(2):
				(*results_file) << boost::lexical_cast<std::string>(r_MacNode_location[0]) << ","
								<< boost::lexical_cast<std::string>(r_MacNode_location[1]) << ",";
				displacement = sqrt(r_MacNode_location[0]*r_MacNode_location[0] + r_MacNode_location[1]*r_MacNode_location[1]);
		break;
		case(3):
				(*results_file) << boost::lexical_cast<std::string>(r_MacNode_location[0]) << ","
								<< boost::lexical_cast<std::string>(r_MacNode_location[1]) << ","
								<< boost::lexical_cast<std::string>(r_MacNode_location[2]) << ",";
				displacement = sqrt(r_MacNode_location[0]*r_MacNode_location[0] + r_MacNode_location[1]*r_MacNode_location[1] + r_MacNode_location[2]*r_MacNode_location[2]);
		break;
		}
		(*results_file) << boost::lexical_cast<std::string>(displacement) << ","
						<< boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTime());
		// Tidy up
		results_file->close();


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
    													}


};


#endif /* TESTSPHEROIDADDINGMACROPHAGES_HPP_ */

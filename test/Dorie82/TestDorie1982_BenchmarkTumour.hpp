#ifndef TESTDORIE1982_BENCHMARKTUMOUR_HPP_
#define TESTDORIE1982_BENCHMARKTUMOUR_HPP_


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
#include "NodeVelocityWriter.hpp"
#include "BoundaryNodeWriter.hpp"

#include "SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic.hpp"
#include "UniformCellCycleModelWithQuiescence.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "ChemotacticForceCSF1.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "EllipticGrowingDomainPdeModifierVariableTimestep.hpp"
#include "RandomNumberGenerator.hpp"
//#include "Debug.hpp"

#include "AveragedSourceEllipticPde.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticBoxDomainPdeModifierVariableTimestep.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ApoptoticCellKiller.hpp"

#include "DiffusionForceChooseD.hpp"
#include "BiasedBrownianMotion.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellVolumesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "UniformCellCycleModelWithQuiescence.hpp"

#include "CuboidPeriodicBoundaryCondition.hpp"
#include "OutputOnlyMacrophageSummaryStatisticsModifier.hpp"
#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis.hpp"
#include "ActualLinearSpringForce.hpp"

#include "RepulsionForce.hpp"
#include "VaryDampingCoefficientOfBoundaryNodes.hpp"
#include "ExternalPressureForce.hpp"
#include "ExternalPressureForceOnSurroundingSphere.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



class TestDorie1982_BenchmarkTumour : public AbstractCellBasedTestSuite
{
public:

	void dontTestSpheroidGrowthCreateBenchmarkTumour() throw(Exception)
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

		double cubeDomainDistanceToBoundary = 30;

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

		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

		// Create tumour cells manually
		for (unsigned i=0; i<nodeNum; i++)
		{
			SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic* p_model = new SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic;
			//TODO Archiving problem
			//			NoCellCycleModel* p_model = new NoCellCycleModel;
			p_model->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1.0);

			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);

			p_model->SetHypoxicConcentration(0.15);
			p_model->SetQuiescentConcentration(0.35);
			p_model->SetCriticalHypoxicDuration(8);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(  p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
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
		int updateIntervalForPdeInTimesteps = 60;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("BenchmarkTumour/120hours");
		simulator.SetSamplingTimestepMultiple(120);
		simulator.SetEndTime(120.0);

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.05);
		simulator.AddForce(p_diffusion_force);

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

		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void dontTestSpheroidGrowthLoadBenchmarkTumour() throw(Exception)
	    																									{
		/* This is a simulation designed to provide groundwork to reproducing Dorie et al.
		 * "Migration and Internalization of Cells and Polystyrene Microspheres in Tumour Cell Spheroids"
		 * Experimental Cell Research 141 (1982) 201-209
		 *
		 * In this case, inert polystyrene beads infiltrate a tumour spheroid where the only forces at work are proliferation and brownian motion
		 */

		OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("BenchmarkTumour/120hours", 240.0);
		p_simulator->SetEndTime(480.0);

		// We remake the cuboid for the PDE FE Mesh
		double cubeDomainDistanceToBoundary = 30;
		ChastePoint<3> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Now loop over simulation modifiers
		for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
				iter != p_simulator->GetSimulationModifiers()->end();
				++iter)
		{

			// Attempt to convert modifier to EllipticBoxDomainPdeModifier. If modifier is of a different type, this should return NULL
			boost::shared_ptr<EllipticBoxDomainPdeModifierVariableTimestep<3> > p_newModifier = boost::dynamic_pointer_cast<EllipticBoxDomainPdeModifierVariableTimestep<3> >(*iter);

			// If not NULL, we have found the correct modifier
			if(p_newModifier)
			{
				// Now regenerate the FeMesh.
				p_newModifier->GenerateFeMesh(p_cuboid,2.0);
				p_newModifier->SetTimestepInterval(60);
			}
		}

		p_simulator->Solve();
		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);
		delete p_simulator;



	    																									}

	void dontTestSpheroidGrowthCreateBenchmarkTumourNecrosisPersists() throw(Exception)
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

		double cubeDomainDistanceToBoundary = 30;

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

		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

		// Create tumour cells manually
		for (unsigned i=0; i<nodeNum; i++)
		{
			SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic* p_model = new SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic;
			p_model->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1.0);
			p_cell->SetApoptosisTime(0.5); // Apoptosis time in hours - how long before a cell is removed?

			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);

			p_model->SetHypoxicConcentration(0.15);
			p_model->SetQuiescentConcentration(0.35);
			p_model->SetCriticalHypoxicDuration(8);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(  p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
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
		int updateIntervalForPdeInTimesteps = 60;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("BenchmarkTumour/NecrosisPersists_0.5HoursOnly");
		simulator.SetSamplingTimestepMultiple(120);
		simulator.SetEndTime(120.0);

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.05);
		simulator.AddForce(p_diffusion_force);

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

		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
									}

	void dontTestSpheroidGrowthLoadBenchmarkTumourNecrosisPersists() throw(Exception)
	    																									{
		/* This is a simulation designed to provide groundwork to reproducing Dorie et al.
		 * "Migration and Internalization of Cells and Polystyrene Microspheres in Tumour Cell Spheroids"
		 * Experimental Cell Research 141 (1982) 201-209
		 *
		 * In this case, inert polystyrene beads infiltrate a tumour spheroid where the only forces at work are proliferation and brownian motion
		 */

		OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("BenchmarkTumour/NecrosisPersists_0.5HoursOnly", 120.0);
		p_simulator->SetEndTime(240.0);

		// We remake the cuboid for the PDE FE Mesh
		double cubeDomainDistanceToBoundary = 30;
		ChastePoint<3> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Now loop over simulation modifiers
		for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
				iter != p_simulator->GetSimulationModifiers()->end();
				++iter)
		{

			// Attempt to convert modifier to EllipticBoxDomainPdeModifier. If modifier is of a different type, this should return NULL
			boost::shared_ptr<EllipticBoxDomainPdeModifierVariableTimestep<3> > p_newModifier = boost::dynamic_pointer_cast<EllipticBoxDomainPdeModifierVariableTimestep<3> >(*iter);

			// If not NULL, we have found the correct modifier
			if(p_newModifier)
			{
				// Now regenerate the FeMesh.
				p_newModifier->GenerateFeMesh(p_cuboid,2.0);
				p_newModifier->SetTimestepInterval(60);
			}
		}

		p_simulator->Solve();
		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);
		delete p_simulator;



	    																									}

	void dontTestSpheroidGrowthCreateBenchmarkTumourVaryDt() throw(Exception)
									{
		/* This is a simulation designed to provide groundwork to reproducing Dorie et al.
		 * "Migration and Internalization of Cells and Polystyrene Microspheres in Tumour Cell Spheroids"
		 * Experimental Cell Research 141 (1982) 201-209
		 *
		 * In this case, inert polystyrene beads infiltrate a tumour spheroid where the only forces at work are proliferation and brownian motion
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;
		const int DIM = 3;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 40;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 1.5;//13.5;
		for (double x=-initialRadius; x<initialRadius; x++)
		{
			for (double y=-initialRadius; y<initialRadius; y++)
			{
				for (double z=-initialRadius; z<initialRadius; z++)
				{
					if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
					{
						nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y, z));
						nodeNum++;
					}
				}
			}
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

		double averageCellCycleLength = 12.0;
		// Create tumour cells manually
		for (unsigned i=0; i<nodeNum; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength-1.0);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength+1.0);
			p_model->SetHypoxicConcentration(0.15);
			p_model->SetQuiescentConcentration(0.35);
			p_model->SetCriticalHypoxicDuration(8);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->SetApoptosisTime(1.0); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength-1.0);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);

		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();

		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary,-cubeDomainDistanceToBoundary);
		ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary,cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		int updateIntervalForPdeInTimesteps = 60;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("BenchmarkTumour/500hours");
		simulator.SetSamplingTimestepMultiple(1);
		simulator.SetEndTime(500.0);
		// Dt default is 30 secs - 0.00833 hours (1.0/120.0)
		double minutesPerTimestep = 6;

		double newTimestep = minutesPerTimestep*(2.0*1.0/120.0);
		simulator.SetDt(newTimestep);
		cell_population.SetAbsoluteMovementThreshold(2.0*2*minutesPerTimestep); // As we've increased the timestep, increase the Movement threshold (default 2.0)


		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<DIM>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		/*// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.05);
		simulator.AddForce(p_diffusion_force);
		 */
		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(3.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();

		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
									}

	void dontTestSpheroidGrowthCreateBenchmarkTumourChangeTimeInterpretation() throw(Exception)
									{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		/* We rescale our time units to effectively increase the timestep
		 */
		int timestepsPerHour = 120; // One timestep is interpreted as every 15 minutes rather than 30 seconds
		int visualisationOutputFrequencyPerHour = 2;


		double chasteDefaultTimestepInHours = 1.0/120.0;
		double intendedNewTimestepInHours = 1.0/timestepsPerHour;
		double timeRescalingConstant = chasteDefaultTimestepInHours/intendedNewTimestepInHours;

		const int DIM = 3;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 40;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 6;
		for (double x=-initialRadius; x<initialRadius+1; x++)
		{
			for (double y=-initialRadius; y<initialRadius+1; y++)
			{
				for (double z=-initialRadius; z<initialRadius+1; z++)
				{
					if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
					{
						nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y, z));
						nodeNum++;
					}
				}
			}
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

		double averageCellCycleLength = 32.0*timeRescalingConstant;
		// Create tumour cells manually
		for (unsigned i=0; i<nodeNum; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*0.9);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*1.1);
			p_model->SetHypoxicConcentration(0.2);
			p_model->SetQuiescentConcentration(0.4);
			p_model->SetCriticalHypoxicDuration(8*timeRescalingConstant);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->SetApoptosisTime(24*timeRescalingConstant); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.9);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(10.0);//Set big movement threshold

		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddPopulationWriter<NodeVelocityWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();

		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary,-cubeDomainDistanceToBoundary);
		ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary,cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		int updateIntervalForPdeInTimesteps = timestepsPerHour/2;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/BenchmarkTumours/NoQuiescentCompartmentDA3";//"BenchmarkTumour/ActualLinearSpringForce";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(timestepsPerHour/visualisationOutputFrequencyPerHour);
		simulator.SetEndTime(200*timeRescalingConstant);

		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(0.3);
		p_macStats_modifier->SetHypoxiaLevel(0.3);
		p_macStats_modifier->SetOutputFrequencyInHours(0.5*timeRescalingConstant); // Every 30 minutes
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.AddSimulationModifier(p_macStats_modifier);

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<DIM>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		//MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		//p_force->SetMeinekeSpringStiffness(5.0);
		//p_force->SetCutOffLength(1.5);
		//simulator.AddForce(p_force);

		//		MAKE_PTR(ActualLinearSpringForce<DIM>, p_force);
		//		p_force->SetMeinekeSpringStiffness(4.0);
		//		p_force->SetCutOffLength(1.5);

		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(4.0);
		p_force->SetCutOffLength(1.5);
		p_force->SetHeterotypicLabelledSpringConstantMultiplier(0.5);
		p_force->SetHomotypicLabelledSpringConstantMultiplier(0.4);
		simulator.AddForce(p_force);


		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();

		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
									}

	void dontTestGrowingDomain() throw(Exception)
									{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		/* We rescale our time units to effectively increase the timestep
		 */
		const int DIM = 2;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		//		double cubeDomainDistanceToBoundary = 40;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 9;
		for (double x=-initialRadius; x<initialRadius+1; x++)
		{

			for (double y=-initialRadius; y<initialRadius+1; y++)
			{

				if(pow(x,2) + pow(y,2) < pow(initialRadius,2))
				{
					nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));
					nodeNum++;
				}

				//											for (double z=-initialRadius; z<initialRadius+1; z++)
				//											{
				//												if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
				//												{
				//													nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y, z));
				//													nodeNum++;
				//												}
				//											}
			}
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

		double averageCellCycleLength = 16;
		// Create tumour cells manually
		double quiescence = 0.4;
		double hypoxia = 0.2;
		for (unsigned i=0; i<nodeNum; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*0.9);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*1.1);
			p_model->SetHypoxicConcentration(hypoxia);
			p_model->SetQuiescentConcentration(quiescence);
			p_model->SetCriticalHypoxicDuration(4);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 0.1);
			p_cell->SetApoptosisTime(24); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.9);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(1000.0);//Set big movement threshold

		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddPopulationWriter<NodeVelocityWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();

		//
		// Oxygen PDE
		//
		MAKE_PTR_ARGS(CellwiseSourceEllipticPde<DIM>, p_pde, (cell_population, -0.1));
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.0));
		bool is_neumann_bc = false;

		MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc));
		p_pde_modifier->SetDependentVariableName("oxygen");
		//		p_pde_modifier->SetTimestepInterval(60); // Update every 60 timesteps = 30 minutes

		//				// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		//				c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
		//				double cubeDomainDistanceToBoundary = 2;
		//				ChastePoint<DIM> lower(-1000, -cubeDomainDistanceToBoundary);
		//				ChastePoint<DIM> upper(1000, cubeDomainDistanceToBoundary);
		//				MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		//				// Create a PDE modifier and set the name of the dependent variable in the PDE
		//				int updateIntervalForPdeInTimesteps = timestepsPerHour/2;
		//				MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifierVariableTimestep<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		//				p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		//				p_pde_modifier->SetDependentVariableName("oxygen");
		//				p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/BenchmarkTumours/GrowingDomain_13_Mar";//"BenchmarkTumour/ActualLinearSpringForce";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(1000);

		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescence);
		p_macStats_modifier->SetHypoxiaLevel(hypoxia);
		p_macStats_modifier->SetOutputFrequencyInHours(0.5); // Every 30 minutes
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.AddSimulationModifier(p_macStats_modifier);

		//				// Add periodic boundary conditions for cells
		//				MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<DIM>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		//				simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);

		//		c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
		//		MAKE_PTR_ARGS(SphereGeometryBoundaryCondition<DIM>,p_sphere_boundary_condition,(&cell_population,centroid,5,5));
		//		simulator.AddCellPopulationBoundaryCondition(p_sphere_boundary_condition);


		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		//MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		//p_force->SetMeinekeSpringStiffness(5.0);
		//p_force->SetCutOffLength(1.5);
		//simulator.AddForce(p_force);

		//		MAKE_PTR(ActualLinearSpringForce<DIM>, p_force);
		//		p_force->SetMeinekeSpringStiffness(4.0);
		//		p_force->SetCutOffLength(1.5);

		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		p_force->SetHeterotypicLabelledSpringConstantMultiplier(0.5);
		p_force->SetHomotypicLabelledSpringConstantMultiplier(0.5);
		p_force->SetCloserThanRestLengthSpringConstantMultiplier(1.1);
		simulator.AddForce(p_force);


		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();

		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
									}


	void dontTestGrowingDomainWithVaryingDampingCoefficient() throw(Exception)
									{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		/* We rescale our time units to effectively increase the timestep
		 */
		const int DIM = 2;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		//		double cubeDomainDistanceToBoundary = 40;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 5;
		for (double x=-initialRadius; x<initialRadius+1; x++)
		{

			for (double y=-initialRadius; y<initialRadius+1; y++)
			{

				if(pow(x,2) + pow(y,2) < pow(initialRadius,2))
				{
					nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));
					nodeNum++;
				}

				//											for (double z=-initialRadius; z<initialRadius+1; z++)
				//											{
				//												if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
				//												{
				//													nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y, z));
				//													nodeNum++;
				//												}
				//											}
			}
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

		double averageCellCycleLength = 24;
		// Create tumour cells manually
		double quiescence = 0.98;
		double hypoxia = 0.95;
		for (unsigned i=0; i<nodeNum; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*0.75);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*1.25);
			p_model->SetHypoxicConcentration(hypoxia);
			p_model->SetQuiescentConcentration(quiescence);
			p_model->SetCriticalHypoxicDuration(0);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1.0);
			p_cell->SetApoptosisTime(12); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.75);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}


		// Make cell population (2D)
		//NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold


		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddPopulationWriter<NodeVelocityWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		cell_population.AddPopulationWriter<BoundaryNodeWriter>();

		//
		// Oxygen PDE
		//
		MAKE_PTR_ARGS(CellwiseSourceEllipticPde<DIM>, p_pde, (cell_population, -0.03));
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.0));
		bool is_neumann_bc = false;

		MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc));
		p_pde_modifier->SetDependentVariableName("oxygen");
		//		p_pde_modifier->SetTimestepInterval(60); // Update every 60 timesteps = 30 minutes


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/BenchmarkTumours/ExternalPressure_26MarHardCoreTEMP";//"BenchmarkTumour/ActualLinearSpringForce";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(60);
		simulator.SetEndTime(1000);

		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescence);
		p_macStats_modifier->SetHypoxiaLevel(hypoxia);
		p_macStats_modifier->SetOutputFrequencyInHours(0.5); // Every 30 minutes
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.AddSimulationModifier(p_macStats_modifier);


		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);


		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(1.0);
		p_force->SetCutOffLength(1.5);
		p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
		p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
		p_force->SetCloserThanRestLengthSpringConstantMultiplier(1.0);
		simulator.AddForce(p_force);


		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);


		//		// Finally, add the "force" which rescales the damping coefficient on the boundary nodes
		//		MAKE_PTR(VaryDampingCoefficientOfBoundaryNodes<DIM>, p_dampingRescale);
		//		p_dampingRescale->SetInteriorNodeDampingConstant(1.0);
		//		p_dampingRescale->SetBoundaryNodeDampingConstant(100.0);
		//		simulator.AddForce(p_dampingRescale);
		MAKE_PTR(ExternalPressureForce<DIM>, p_pressure);
		p_pressure->SetPressure(1.0);
		simulator.AddForce(p_pressure);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();

		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
									}

	void TestSpheroidGrowthBoxWithInwardForce() throw(Exception)
										{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		const int DIM = 2;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 40;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 5;
		for (double x=-initialRadius; x<initialRadius+1; x++)
		{

			for (double y=-initialRadius; y<initialRadius+1; y++)
			{

				if(pow(x,2) + pow(y,2) < pow(initialRadius,2))
				{
					nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));
					nodeNum++;
				}

//															for (double z=-initialRadius; z<initialRadius+1; z++)
//															{
//																if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
//																{
//																	nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y, z));
//																	nodeNum++;
//																}
//															}
			}
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

		double averageCellCycleLength = 24;
		// Create tumour cells manually
		double quiescence = 0.5;
		double hypoxia = 0.4;
		for (unsigned i=0; i<nodeNum; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*0.75);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*1.25);
			p_model->SetHypoxicConcentration(hypoxia);
			p_model->SetQuiescentConcentration(quiescence);
			p_model->SetCriticalHypoxicDuration(8);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1.0);
			p_cell->SetApoptosisTime(24); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.75);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold


		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddPopulationWriter<NodeVelocityWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		cell_population.AddPopulationWriter<BoundaryNodeWriter>();

		// Make PDE (Oxygen)
		double consumptionRate = -0.02;
		double diffusionCoefficient = 1;
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_pde, (cell_population, consumptionRate,diffusionCoefficient));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary,-cubeDomainDistanceToBoundary);
		ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary,cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		int updateIntervalForPdeInTimesteps = 120/2;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(false); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/BenchmarkTumours/28Mar_IncreasedDiffusionCoefficient2";//"BenchmarkTumour/ActualLinearSpringForce";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(60);
		simulator.SetEndTime(1000);

		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescence);
		p_macStats_modifier->SetHypoxiaLevel(hypoxia);
		p_macStats_modifier->SetOutputFrequencyInHours(0.5); // Every 30 minutes
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.AddSimulationModifier(p_macStats_modifier);

		//			// Add periodic boundary conditions for cells
		//			MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<DIM>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		//			simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		//MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		//p_force->SetMeinekeSpringStiffness(5.0);
		//p_force->SetCutOffLength(1.5);
		//simulator.AddForce(p_force);

		//		MAKE_PTR(ActualLinearSpringForce<DIM>, p_force);
		//		p_force->SetMeinekeSpringStiffness(4.0);
		//		p_force->SetCutOffLength(1.5);

		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(4.0);
		p_force->SetCutOffLength(1.5);
		p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
		p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
		simulator.AddForce(p_force);


		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		//		// Finally, add the "force" which rescales the damping coefficient on the boundary nodes
		//		MAKE_PTR(VaryDampingCoefficientOfBoundaryNodes<DIM>, p_dampingRescale);
		//		p_dampingRescale->SetInteriorNodeDampingConstant(1.0);
		//		p_dampingRescale->SetBoundaryNodeDampingConstant(100.0);
		//		simulator.AddForce(p_dampingRescale);
		//MAKE_PTR(ExternalPressureForce<DIM>, p_pressure);
		MAKE_PTR(ExternalPressureForceOnSurroundingSphere<DIM>, p_pressure);
		p_pressure->SetPressure(10);
		simulator.AddForce(p_pressure);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();

		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
										}


};


#endif /* TESTDORIE1982_BENCHMARKTUMOUR_HPP_ */

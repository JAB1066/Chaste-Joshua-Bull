#ifndef TESTDORIE1982_ADDINGMACROPHAGES_HPP_
#define TESTDORIE1982_ADDINGMACROPHAGES_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"

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

#include "AveragedSourceEllipticPde.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "AveragedSourceEllipticPdeOxygenBelowThreshold.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticBoxDomainPdeModifierVariableTimestep.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ApoptoticCellKiller.hpp"

#include "DiffusionForceChooseD.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

//#include "NodeLocationWriter.hpp"
#include "MacrophageNodeLocationWriter.hpp"
#include "UniformCellCycleModelWithQuiescence.hpp"

#include "CuboidPeriodicBoundaryCondition.hpp"

#include "AddMacrophagesAtSpecifiedSpheroidSizeModifier.hpp"
#include "AddMacrophagesAtSpecifiedTimeModifier.hpp"
#include "AddMacrophagesToBoundaryNodesAtSpecifiedTimeModifier.hpp"

#include "ExternalPressureForceOnSurroundingSphere.hpp"
#include "ExternalPressureForceOnConcaveHull.hpp"

#include "OutputOnlyMacrophageSummaryStatisticsModifier.hpp"
#include "OutputOnlyMacrophageSummaryStatisticsModifierWithCSF1.hpp"
#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis.hpp"
#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages.hpp"
#include "NodeLocationWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "BoundaryNodeWriter.hpp"
#include "BoundaryCellWriter.hpp"

#include "ChemotacticForceCSF1.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>
#include "Debug.hpp"


class TestDorie1982_AddingMacrophages : public AbstractCellBasedTestSuite
{
public:

	void dontTestMacrophageCellWriter() throw(Exception)
	{
		/*
		 * In this test, "macrophages" have no macrophage properties - they are completely inert,
		 * with the exception of when a degree of Brownian Motion applies to the entire simulation.
		 *
		 * As a result, we can interpret them as inert "microspheres" a la Dorie 1982
		 *
		 */

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 30;

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);


		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 3;//13.5;
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

			double hypoxicConcentration = 0.15;
			double quiescentConcentration = 0.35;
			p_model->SetHypoxicConcentration(hypoxicConcentration);
			p_model->SetQuiescentConcentration(quiescentConcentration);
			p_model->SetCriticalHypoxicDuration(8);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(  p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}



		/*
		 * Now we add our macrophages
		 */
		unsigned firstMacrophageNode = nodes.size();

		// Macrophage Nodes - in a shell at fixed radius
		// Number of Macrophages remains constant throughout simulation
		unsigned numMacrophages=20;

		std::vector<double> xVec;
		std::vector<double> yVec;
		std::vector<double> zVec;
		unsigned macCount=0;
		double random_x;
		double random_y;
		double zCoord;

		// Add macrophages at random points at edge of spheroid
		double macrophageSphereRadius = 4;
		while(macCount < numMacrophages)
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			random_x = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			random_y = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));

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
			macCount++;
			nodes.push_back(new Node<3>(nodeNum,  false,  random_x, random_y, zCoord));
			nodeNum++;
		}


		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


		// Now create macrophages
		for (unsigned i=firstMacrophageNode; i<nodes.size(); i++)
		{
			// Make Macrophage
			NoCellCycleModel* p_model2 = new NoCellCycleModel;
			p_model2->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model2));
			p_cell->SetCellProliferativeType(p_macrophage_type);

			cells.push_back(p_cell);
		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		//cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddPopulationWriter<MacrophageNodeLocationWriter>();

		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
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
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
		int updateIntervalForPdeInTimesteps = 60;
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/AddingMacrophages/CellWriterTest";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(60); // One visualisation every 30 minutes...
		double simDuration = 2.0;
		simulator.SetEndTime(simDuration); // ...for 2 hours

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.02);
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

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void dontTestMacrophageSummaryStatsModifier() throw(Exception)
															{
		/*
		 * In this test, "macrophages" have no macrophage properties - they are completely inert,
		 * with the exception of when a degree of Brownian Motion applies to the entire simulation.
		 *
		 * As a result, we can interpret them as inert "microspheres" a la Dorie 1982
		 *
		 */

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 30;
		double hypoxicConcentration = 0.15;
		double quiescentConcentration = 0.35;

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);


		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 3;//13.5;
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
			p_model->SetHypoxicConcentration(hypoxicConcentration);
			p_model->SetQuiescentConcentration(quiescentConcentration);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(  p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}



		/*
		 * Now we add our macrophages
		 */
		unsigned firstMacrophageNode = nodes.size();

		// Macrophage Nodes - in a shell at fixed radius
		// Number of Macrophages remains constant throughout simulation
		unsigned numMacrophages=20;

		std::vector<double> xVec;
		std::vector<double> yVec;
		std::vector<double> zVec;
		unsigned macCount=0;
		double random_x;
		double random_y;
		double zCoord;

		// Add macrophages at random points at edge of spheroid
		double macrophageSphereRadius = 4;
		while(macCount < numMacrophages)
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			random_x = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			random_y = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));

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
			macCount++;
			nodes.push_back(new Node<3>(nodeNum,  false,  random_x, random_y, zCoord));
			nodeNum++;
		}
		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


		// Now create macrophages
		for (unsigned i=firstMacrophageNode; i<nodes.size(); i++)
		{
			// Make Macrophage
			NoCellCycleModel* p_model2 = new NoCellCycleModel;
			p_model2->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model2));
			p_cell->SetCellProliferativeType(p_macrophage_type);

			cells.push_back(p_cell);
		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		//cell_population.AddPopulationWriter<NodeLocationWriter>();
		//cell_population.AddPopulationWriter<MacrophageNodeLocationWriter>();

		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
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
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
		int updateIntervalForPdeInTimesteps = 60;
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<3>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescentConcentration);
		p_macStats_modifier->SetHypoxiaLevel(hypoxicConcentration);
		p_macStats_modifier->SetOutputFrequencyInHours(0.1);

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.AddSimulationModifier(p_macStats_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/AddingMacrophages/StatsModifierTest";
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(UINT_MAX); // One visualisation every 30 minutes...
		double simDuration = 2.0;
		simulator.SetEndTime(simDuration); // ...for 2 hours

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.02);
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
		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
															}

	void dontTestMacrophageSummaryStatsModifierBigSpheroid() throw(Exception)
															{
		/*
		 * In this test, "macrophages" have no macrophage properties - they are completely inert,
		 * with the exception of when a degree of Brownian Motion applies to the entire simulation.
		 *
		 * As a result, we can interpret them as inert "microspheres" a la Dorie 1982
		 *
		 */

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 30;
		double hypoxicConcentration = 0.4;
		double quiescentConcentration = 0.5;

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		/*
		 * We load in another simulation so that we can start with a "burnt in" tumour spheroid
		 */
		OffLatticeSimulation<3>* p_archivedSimulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("BenchmarkTumour/120hours", 120.0);

		NodeBasedCellPopulation<3>& r_archivedPopulation = dynamic_cast<NodeBasedCellPopulation<3> &>(p_archivedSimulator->rGetCellPopulation());

		std::list<CellPtr> archived_cells = r_archivedPopulation.rGetCells();
		std::list<CellPtr>::iterator it;
		int nodeNum = 0;
		for (it = archived_cells.begin(); it != archived_cells.end(); it++)
		{
			// Give each cell a pointer to the property registry (we have taken ownership in this constructor)
			c_vector<double, 3> location = r_archivedPopulation.GetLocationOfCellCentre(*it);
			nodes.push_back(new Node<3>(nodeNum,  false,  location[0], location[1], location[2]));
			double o2conc = (*it)->GetCellData()->GetItem("oxygen");
			double birthTime = (*it)->GetBirthTime();


			SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic* p_model = new SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic;

			p_model->SetDimension(3);
			//p_model->SetStemCellG1Duration(8);
			//p_model->SetTransitCellG1Duration(8);
			p_model->SetHypoxicConcentration(hypoxicConcentration);
			p_model->SetQuiescentConcentration(quiescentConcentration);
			p_model->SetCriticalHypoxicDuration(8);

			CellPtr p_cell(new Cell(p_state, p_model));

			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", o2conc);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

			nodeNum++;
		}


		/*
		 * Now we add our macrophages
		 */
		unsigned firstMacrophageNode = nodes.size();

		// Macrophage Nodes - in a shell at fixed radius
		// Number of Macrophages remains constant throughout simulation
		unsigned numMacrophages=100;

		std::vector<double> xVec;
		std::vector<double> yVec;
		std::vector<double> zVec;
		unsigned macCount=0;
		double random_x;
		double random_y;
		double zCoord;

		// Add macrophages at random points at edge of spheroid
		double macrophageSphereRadius = 9.5;
		while(macCount < numMacrophages)
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			random_x = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
			random_y = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));

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
			macCount++;
			nodes.push_back(new Node<3>(nodeNum,  false,  random_x, random_y, zCoord));
			nodeNum++;
		}


		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


		// Now create macrophages
		for (unsigned i=firstMacrophageNode; i<nodes.size(); i++)
		{
			// Make Macrophage
			NoCellCycleModel* p_model2 = new NoCellCycleModel;
			p_model2->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model2));
			p_cell->SetCellProliferativeType(p_macrophage_type);

			cells.push_back(p_cell);
		}

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(DBL_MAX);
		//cell_population.AddPopulationWriter<NodeLocationWriter>();
		//cell_population.AddPopulationWriter<MacrophageNodeLocationWriter>();

		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
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
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
		int updateIntervalForPdeInTimesteps = 60;
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<3>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescentConcentration);
		p_macStats_modifier->SetHypoxiaLevel(hypoxicConcentration);
		p_macStats_modifier->SetOutputFrequencyInHours(0.1);

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.AddSimulationModifier(p_macStats_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/AddingMacrophages/StatsModifierTest_BigTumour";
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(60);//UINT_MAX); // One visualisation every 30 minutes...
		double simDuration = 96.0;
		simulator.SetEndTime(120.0+simDuration); // ...for 2 hours

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.02);
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
		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
															}

	void dontTestFinalSimulations() throw(Exception)
															{

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		/* We rescale our time units to effectively increase the timestep
		 */
		int timestepsPerHour = 120; // One timestep is interpreted as every 15 minutes rather than 30 seconds
		int visualisationOutputFrequencyPerHour = 2;

		unsigned numMacrophages = 1000;


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
		double initialRadius = 8.0;
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


		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		double averageCellCycleLength = 24.0*timeRescalingConstant;
		// Create tumour cells manually
		for (unsigned i=0; i<nodeNum; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*0.9);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*1.1);
			p_model->SetHypoxicConcentration(0.2);
			p_model->SetQuiescentConcentration(0.4);
			p_model->SetCriticalHypoxicDuration(4*timeRescalingConstant);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->SetApoptosisTime(48*timeRescalingConstant); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.9);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}

		/*
		 * Now we add our macrophages
		 */
		unsigned firstMacrophageNode = nodes.size();

		// Macrophage Nodes - in a shell at fixed radius
		// Number of Macrophages remains constant throughout simulation

		unsigned macCount=0;
		double random_x;
		double random_y;
		double random_z;
		double normalize;

		// Add macrophages at random points at edge of spheroid
		double macrophageSphereRadius = 9.0;
		while(macCount < numMacrophages)
		{
			random_x = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
			random_y = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
			random_z = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();

			normalize = macrophageSphereRadius/sqrt(pow(random_x,2) + pow(random_y,2)+ pow(random_z,2));
			// Make vectors lie on sphere of radiusmacrophageSphereRadius

			nodes.push_back(new Node<3>(nodeNum,  false,  random_x*normalize, random_y*normalize, random_z*normalize));
			macCount++;
			nodeNum++;
		}


		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


		// Now create macrophages
		for (unsigned i=firstMacrophageNode; i<nodes.size(); i++)
		{
			// Make Macrophage
			NoCellCycleModel* p_model2 = new NoCellCycleModel;
			p_model2->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model2));
			p_cell->SetCellProliferativeType(p_macrophage_type);

			cells.push_back(p_cell);
		}


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(100.0);//Set big movement threshold

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
		output_directory << "AddingMacrophages/LPLD_WithVelocityWriter_19Dec";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(timestepsPerHour/visualisationOutputFrequencyPerHour);
		simulator.SetEndTime(240*timeRescalingConstant);

		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(0.4);
		p_macStats_modifier->SetHypoxiaLevel(0.2);
		p_macStats_modifier->SetOutputFrequencyInHours(1*timeRescalingConstant); // Every 30 minutes
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

		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(4.0);
		p_force->SetCutOffLength(1.5);
		//p_force->SetHeterotypicLabelledSpringConstantMultiplier(0.5);
		//p_force->SetHomotypicLabelledSpringConstantMultiplier(0.75);
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

	void dontTestFinalSimulationWithLoadedSpheroid() throw(Exception)
																			{

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		/* We rescale our time units to effectively increase the timestep
		 */
		int timestepsPerHour = 120; // One timestep is interpreted as every 15 minutes rather than 30 seconds
		int visualisationOutputFrequencyPerHour = 2;

		unsigned numMacrophages = 100;


		double chasteDefaultTimestepInHours = 1.0/120.0;
		double intendedNewTimestepInHours = 1.0/timestepsPerHour;
		double timeRescalingConstant = chasteDefaultTimestepInHours/intendedNewTimestepInHours;

		const int DIM = 3;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 40;

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);
		/*
		 * We load in another simulation so that we can start with a "burnt in" tumour spheroid
		 */
		double simulationLoadTime = 120.0;
		OffLatticeSimulation<3>* p_archivedSimulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("BenchmarkTumour/120hours", simulationLoadTime);

		NodeBasedCellPopulation<3>& r_archivedPopulation = dynamic_cast<NodeBasedCellPopulation<3> &>(p_archivedSimulator->rGetCellPopulation());

		std::list<CellPtr> archived_cells = r_archivedPopulation.rGetCells();
		std::list<CellPtr>::iterator it;

		double averageCellCycleLength = 16;
		double hypoxicConcentration = 0.3;
		double quiescentConcentration = 0.3;
		double criticalHypoxicDuration = 16;

		int nodeNum = 0;
		double maxRadius = 0;
		for (it = archived_cells.begin(); it != archived_cells.end(); it++)
		{
			// Give each cell a pointer to the property registry (we have taken ownership in this constructor)
			c_vector<double, 3> location = r_archivedPopulation.GetLocationOfCellCentre(*it);
			nodes.push_back(new Node<3>(nodeNum,  false,  location[0], location[1], location[2]));
			if(pow(location[0],2) + pow(location[1],2) + pow(location[2],2) > maxRadius)
			{
				maxRadius = pow(location[0],2) + pow(location[1],2) + pow(location[2],2);
			}

			double o2conc = (*it)->GetCellData()->GetItem("oxygen");
			double birthTime = simulationLoadTime - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.9);

			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(3);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*0.9);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*1.1);
			p_model->SetHypoxicConcentration(hypoxicConcentration);
			p_model->SetQuiescentConcentration(quiescentConcentration);
			p_model->SetCriticalHypoxicDuration(criticalHypoxicDuration);
			p_model->SetBirthTime(birthTime);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", o2conc);
			p_cell->SetApoptosisTime(48.0); // Apoptosis time in hours - how long before a cell is removed?
			p_cell->SetCellCycleModel(p_model);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

			nodeNum++;
		}



		/*
		 * Now we add our macrophages
		 */
		unsigned firstMacrophageNode = nodes.size();

		// Macrophage Nodes - in a shell at fixed radius
		// Number of Macrophages remains constant throughout simulation

		unsigned macCount=0;
		double random_x;
		double random_y;
		double random_z;
		double normalize;
		double macrophageSphereRadius =  sqrt(abs(maxRadius))+1;

		// Add macrophages at random points at edge of spheroid
		while(macCount < numMacrophages)
		{
			random_x = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
			random_y = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
			random_z = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();

			normalize = macrophageSphereRadius/sqrt(pow(random_x,2) + pow(random_y,2)+ pow(random_z,2));
			// Make vectors lie on sphere of radiusmacrophageSphereRadius

			nodes.push_back(new Node<3>(nodeNum,  false,  random_x*normalize, random_y*normalize, random_z*normalize));
			macCount++;
			nodeNum++;
		}


		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


		// Now create macrophages
		for (unsigned i=firstMacrophageNode; i<nodes.size(); i++)
		{
			// Make Macrophage
			NoCellCycleModel* p_model2 = new NoCellCycleModel;
			p_model2->SetDimension(3);

			CellPtr p_cell(new Cell(p_state, p_model2));
			p_cell->SetCellProliferativeType(p_macrophage_type);

			cells.push_back(p_cell);
		}


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(100.0);//Set big movement threshold

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
		output_directory << "Dorie1982/ExampleSimulations/NoQuiescentCompartment";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(timestepsPerHour/visualisationOutputFrequencyPerHour);
		simulator.SetEndTime(300*timeRescalingConstant);

		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescentConcentration);
		p_macStats_modifier->SetHypoxiaLevel(hypoxicConcentration);
		p_macStats_modifier->SetOutputFrequencyInHours(1*timeRescalingConstant); // Every 30 minutes
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

		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(4.0);
		p_force->SetCutOffLength(1.5);
		//p_force->SetHeterotypicLabelledSpringConstantMultiplier(0.5);
		//p_force->SetHomotypicLabelledSpringConstantMultiplier(0.75);
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

	void dontTestAddingMacrophagesWhenTumourReachesSpecifiedSize() throw(Exception)
															{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		int visualisationOutputFrequencyPerHour = 2;

		const int DIM = 3;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 40;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 3.0;
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


		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		double averageCellCycleLength = 16.0;
		// Create tumour cells manually
		for (unsigned i=0; i<nodeNum; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*0.9);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*1.1);
			p_model->SetHypoxicConcentration(0.2);
			p_model->SetQuiescentConcentration(0.4);
			p_model->SetCriticalHypoxicDuration(4);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.9);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}



		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(6.0);//Set big movement threshold

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
		int updateIntervalForPdeInTimesteps = 120/2;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		std::stringstream output_directory;
		output_directory << "AddingMacrophages/AddMacrophagesAtSpecifiedSize";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
		simulator.SetEndTime(50);

		// Add macrophages at set time
		MAKE_PTR(AddMacrophagesAtSpecifiedSpheroidSizeModifier<DIM>, p_addMacs_modifier);
		p_addMacs_modifier->SetNumberOfMacrophagesToAdd(10);
		p_addMacs_modifier->SetNumberOfTumourCellsBeforeMacrophagesAreAdded(100);
		p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(10);
		simulator.AddSimulationModifier(p_addMacs_modifier);


		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(0.4);
		p_macStats_modifier->SetHypoxiaLevel(0.2);
		p_macStats_modifier->SetOutputFrequencyInHours(0.5); // Every 30 minutes
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

	void dontTestAddingMacrophagesWhenSimulationReachesSpecifiedTime() throw(Exception)
															{

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		int visualisationOutputFrequencyPerHour = 2;
		const int DIM = 2;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 50;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 12.0;
		for (double x=-initialRadius; x<initialRadius+1; x++)
		{
			for (double y=-initialRadius; y<initialRadius+1; y++)
			{

				if(pow(x,2) + pow(y,2) < pow(initialRadius,2))
				{
					nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));
					nodeNum++;
				}
				//								for (double z=-initialRadius; z<initialRadius+1; z++)
				//								{
				//									if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
				//									{
				//										nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y, z));
				//										nodeNum++;
				//									}
				//								}
			}
		}


		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		double averageCellCycleLength = 24.0;
		// Create tumour cells manually
		double quiescence = 0.7;
		double hypoxia = 0.6;
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
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.75);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}



		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddPopulationWriter<NodeVelocityWriter>();
		cell_population.AddPopulationWriter<BoundaryNodeWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		cell_population.AddCellWriter<BoundaryCellWriter>();
		//cell_population.AddPopulationWriter<BoundaryNodeWriter>();

		// Make PDE (Oxygen)
		double consumptionRate = -0.03;//was -0.03
		double diffusionCoefficient = 1;
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_pde, (cell_population, consumptionRate,diffusionCoefficient));
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.1));
		bool is_neumann_bc = false; // Dirichlet BCs


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary,-cubeDomainDistanceToBoundary);
		ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary,cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		int updateIntervalForPdeInTimesteps = 120/2;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 1.0));
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(false); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/AddingMacrophages/TestBeadVolumeFix";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
		simulator.SetEndTime(400);

		// Add macrophages at set time
		MAKE_PTR(AddMacrophagesAtSpecifiedTimeModifier<DIM>, p_addMacs_modifier);
		p_addMacs_modifier->SetNumberOfMacrophagesToAdd(50);
		p_addMacs_modifier->SetTimeToAddMacrophages(75);
		//p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(250);
		simulator.AddSimulationModifier(p_addMacs_modifier);


		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescence);
		p_macStats_modifier->SetHypoxiaLevel(hypoxia);
		p_macStats_modifier->SetOutputFrequencyInHours(0.5); // Every 30 minutes
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.AddSimulationModifier(p_macStats_modifier);

		//		// Add periodic boundary conditions for cells
		//		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<DIM>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		//		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		//MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		//p_force->SetMeinekeSpringStiffness(5.0);
		//p_force->SetCutOffLength(1.5);
		//simulator.AddForce(p_force);

		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(15.0);
		p_force->SetCutOffLength(1.5);
		p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
		p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
		simulator.AddForce(p_force);


		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		MAKE_PTR(ExternalPressureForceOnConcaveHull<DIM>, p_pressure);
		p_pressure->SetPressure(1);
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

	void TestAddingMacrophagesWithCSF1() throw(Exception)
															{

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		int visualisationOutputFrequencyPerHour = 2;
		const int DIM = 2;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 50;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 6.0;
		for (double x=-initialRadius; x<initialRadius+1; x++)
		{
			for (double y=-initialRadius; y<initialRadius+1; y++)
			{

				if(pow(x,2) + pow(y,2) < pow(initialRadius,2))
				{
					nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));
					nodeNum++;
				}
				//								for (double z=-initialRadius; z<initialRadius+1; z++)
				//								{
				//									if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
				//									{
				//										nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y, z));
				//										nodeNum++;
				//									}
				//								}
			}
		}


		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		double averageCellCycleLength = 16.0;
		// Create tumour cells manually
		double quiescence = 0.7;
		double hypoxia = 0.5;
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
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->GetCellData()->SetItem("csf1", 1.0);
			p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.75);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}



		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddPopulationWriter<NodeVelocityWriter>();
		cell_population.AddPopulationWriter<BoundaryNodeWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		cell_population.AddCellWriter<BoundaryCellWriter>();
		//cell_population.AddPopulationWriter<BoundaryNodeWriter>();

		// Make PDE (Oxygen)
		double consumptionRate = -0.03;//was -0.03
		double diffusionCoefficient = 1;
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_pde, (cell_population, consumptionRate,diffusionCoefficient));
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
		ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary,-cubeDomainDistanceToBoundary);
		ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary,cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		int updateIntervalForPdeInTimesteps = 120/2;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 1.0));
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(false); //was false



		// Make PDE (CSF1)
		// Assume CSF1 diffuses 4 times slower than o2
		MAKE_PTR_ARGS(AveragedSourceEllipticPdeOxygenBelowThreshold<DIM>, p_pdeCSF, (cell_population, 0.01,0.25,hypoxia));
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bcCSF, (0.0));

		bool is_neumann_bcCSF = false; // Dirichlet BCs
		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboidCSF, (lower, upper));
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<DIM>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboidCSF, 1.0));
		p_pde_modifierCSF->SetDependentVariableName("csf1");
		p_pde_modifierCSF->SetOutputGradient(true);
		p_pde_modifierCSF->SetBcsOnBoxBoundary(true); //was false


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);

		simulator.AddSimulationModifier(p_pde_modifierCSF);

		std::stringstream output_directory;
		output_directory << "Dorie1982/AddingMacrophages/TestCSF1_HypoxicOnly";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
		simulator.SetEndTime(400);

		// Add macrophages at set time
		MAKE_PTR(AddMacrophagesAtSpecifiedTimeModifier<DIM>, p_addMacs_modifier);
		p_addMacs_modifier->SetNumberOfMacrophagesToAdd(50);
		p_addMacs_modifier->SetTimeToAddMacrophages(100);
		//p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(250);
		simulator.AddSimulationModifier(p_addMacs_modifier);


		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifierWithCSF1<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescence);
		p_macStats_modifier->SetHypoxiaLevel(hypoxia);
		p_macStats_modifier->SetOutputFrequencyInHours(0.5); // Every 30 minutes
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.AddSimulationModifier(p_macStats_modifier);

		//		// Add periodic boundary conditions for cells
		//		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<DIM>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		//		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);

		// Add Chemotaxis
		MAKE_PTR(ChemotacticForceCSF1<DIM>, p_chemotactic_force);
		simulator.AddForce(p_chemotactic_force);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		//MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		//p_force->SetMeinekeSpringStiffness(5.0);
		//p_force->SetCutOffLength(1.5);
		//simulator.AddForce(p_force);

		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(15.0);
		p_force->SetCutOffLength(1.5);
		p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
		p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
		simulator.AddForce(p_force);


		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		MAKE_PTR(ExternalPressureForceOnConcaveHull<DIM>, p_pressure);
		p_pressure->SetPressure(5);
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


#endif /* TESTDORIE1982_ADDINGMACROPHAGES_HPP_ */

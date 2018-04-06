#ifndef TEST_IMO17_FIBROBLASTTUMOURPATTERNS_HPP_
#define TEST_IMO17_FIBROBLASTTUMOURPATTERNS_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
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

//#include "SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
//#include "ChemotacticForceCSF1.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
//#include "Debug.hpp"

#include "AveragedSourceParabolicPde.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ApoptoticCellKiller.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "CuboidPeriodicBoundaryCondition.hpp"
#include "VolumeTrackingModifier.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "RandomCellKiller.hpp"
#include "AbstractCellProperty.hpp"

#include "MacrophageProximityLabellingModifier.hpp"
#include "AddProliferationInhibitingTreatmentModifier.hpp"
#include "UniformContactInhibitionCellCycleModelSlowByFibroblasts.hpp"
#include "ResistantIMOMutationState.hpp"
#include "FibroblastIMOMutationState.hpp"
//#include "ResistantMutationStateTrackingModifier.hpp"
#include "RandomCellKillerIgnoringFibroblastsIMO.hpp"
#include "KeepMacrophagesStationaryLinearSpringForce.hpp"
#include "RepulsionForceKeepMacrophagesStationary.hpp"

#include "ColumnDataReader.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



class TestIMO17_FibroblastTumourPatterns : public AbstractCellBasedTestSuite
{
public:

	void dontTestContactInhibition() throw(Exception)
	{
		/* We throw down fibroblasts randomly across the domain, and seed
		 * it with a tumour cell population. Fibroblasts produce growth factor GF, which
		 * promotes the growth of tumour cells.
		 *
		 * After time T, we simulate treatment: a "resistant" TC population is unaffected.
		 * The "non-resistant" cell population only proliferates in the presence of GF.
		 *
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Parameters
		unsigned numberOfFibroblasts = 50;
		unsigned numberOfTumourCells = 500;
		double squareDomainWidth = 20;
		const int DIM = 2;

		//double squareDomainWidthHalved = squareDomainWidth*0.5;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		// We randomly distribute nodes across the domain
		unsigned nodeNum=0;
		double x;
		double y;
		//double z;

		while(nodeNum < (numberOfFibroblasts + numberOfTumourCells))
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			x = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			y = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			//z = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf()-squareDomainWidthHalved);

			nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));//, z));
			nodeNum++;
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_fibroblast_type);

		// Create tumour cells manually
		for (unsigned i=0; i<numberOfTumourCells; i++)
		{
			ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetTransitCellG1Duration(1.0);
			p_model->SetStemCellG1Duration(1.0);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					( p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
		for (unsigned i = 0; i<numberOfFibroblasts; i++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_fibroblast_type);

            cells.push_back(p_cell);
        }


		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);

		// Make PDE for HGF
		double duDtCoefficient=1.0;
        double diffusionCoefficient=1.0;
        double sourceCoefficient=0.03;
		MAKE_PTR_ARGS(CellwiseSourceParabolicPde<DIM>,p_pde,(cell_population, duDtCoefficient,diffusionCoefficient,sourceCoefficient));
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.0));
		//bool is_neumann_bc = true;

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		ChastePoint<DIM> lower(0,0);
		ChastePoint<DIM> upper(squareDomainWidth,squareDomainWidth);
		//ChastePoint<3> lower(-squareDomainWidthHalved, -squareDomainWidthHalved);//, -squareDomainWidthHalved-0.1);
		//ChastePoint<3> upper(squareDomainWidthHalved, squareDomainWidthHalved);//, squareDomainWidthHalved+0.1);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		/*
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifierVariable<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("HGF");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
        */

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		//simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("IMO2017/ContactInhibition/");
		simulator.SetSamplingTimestepMultiple(12); // 12 timesteps = 6 minutes = 0.1 hours
		simulator.SetEndTime(60.0);

		// Add boundary conditions
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = unit_vector<double>(2,1);
        normal(0) = -1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc1, (&cell_population, point, normal));
        p_bc1->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        normal(0) = 0.0; normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc2, (&cell_population, point, normal));
        p_bc2->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = squareDomainWidth; point(1) = squareDomainWidth;
        normal(0) = 0.0; normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc3, (&cell_population, point, normal));
        p_bc3->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        normal(0) = 1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc4, (&cell_population, point, normal));
        p_bc4->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Add a cell killer
        double probabiltyOfDeathInOneHour = 0.01;
        MAKE_PTR_ARGS(RandomCellKiller<DIM>, p_random_killer, (&cell_population, probabiltyOfDeathInOneHour));
        simulator.AddCellKiller(p_random_killer);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(15.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();

		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void dontTestAddingDifferentialProliferationRateByFibroblasts() throw(Exception)
	{
		/* We throw down fibroblasts randomly across the domain, and seed
		 * it with a tumour cell population. Fibroblasts produce growth factor GF, which
		 * promotes the growth of tumour cells.
		 *
		 * After time T, we simulate treatment: a "resistant" TC population is unaffected.
		 * The "non-resistant" cell population only proliferates in the presence of GF.
		 *
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Parameters
		unsigned numberOfFibroblasts = 5;
		unsigned numberOfTumourCells = 500;
		double squareDomainWidth = 20;
		const int DIM = 2;

		//double squareDomainWidthHalved = squareDomainWidth*0.5;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		// We randomly distribute nodes across the domain
		unsigned nodeNum=0;
		double x;
		double y;
		//double z;

		while(nodeNum < (numberOfFibroblasts + numberOfTumourCells))
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			x = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			y = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			//z = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf()-squareDomainWidthHalved);

			nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));//, z));
			nodeNum++;
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_fibroblast_type);

		// Create tumour cells manually
		for (unsigned i=0; i<numberOfTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(0.1);
			p_model->SetMaxCellCycleDuration(1.0);
			p_model->SetMinCellCycleDurationByFibroblast(0.1);
			p_model->SetMaxCellCycleDurationByFibroblast(1.0);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					( p_model->GetMinCellCycleDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
		for (unsigned i = 0; i<numberOfFibroblasts; i++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_fibroblast_type);

            cells.push_back(p_cell);
        }


		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		ChastePoint<DIM> lower(0,0);
		ChastePoint<DIM> upper(squareDomainWidth,squareDomainWidth);
		//ChastePoint<3> lower(-squareDomainWidthHalved, -squareDomainWidthHalved);//, -squareDomainWidthHalved-0.1);
		//ChastePoint<3> upper(squareDomainWidthHalved, squareDomainWidthHalved);//, squareDomainWidthHalved+0.1);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		/*
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifierVariable<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("HGF");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
        */

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		//simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("IMO2017/AddingDifferentialProliferation/");
		simulator.SetSamplingTimestepMultiple(12); // 12 timesteps = 6 minutes = 0.1 hours
		simulator.SetEndTime(24.0);



		// Add boundary conditions
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = unit_vector<double>(2,1);
        normal(0) = -1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc1, (&cell_population, point, normal));
        p_bc1->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        normal(0) = 0.0; normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc2, (&cell_population, point, normal));
        p_bc2->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = squareDomainWidth; point(1) = squareDomainWidth;
        normal(0) = 0.0; normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc3, (&cell_population, point, normal));
        p_bc3->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        normal(0) = 1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc4, (&cell_population, point, normal));
        p_bc4->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Add a cell killer
        double probabiltyOfDeathInOneHour = 0.01*60;
        MAKE_PTR_ARGS(RandomCellKiller<DIM>, p_random_killer, (&cell_population, probabiltyOfDeathInOneHour));
        simulator.AddCellKiller(p_random_killer);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(15.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Label cells which are in close proximity to fibroblasts
		MAKE_PTR(MacrophageProximityLabellingModifier<DIM>, p_proximity_labeller);
		p_proximity_labeller->SetMacrophageProximityLabellingRadius(3.0);
		simulator.AddSimulationModifier(p_proximity_labeller);


		/*** To run the simulation, we call {{{Solve()}}}. ***/
		simulator.Solve();

		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

		/*** We now need to add a simulation modifier to simulate adding the treatment */



		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void dontTestAddingResistantTumourCellPopulation() throw(Exception)
	{
		/* We throw down fibroblasts randomly across the domain, and seed
		 * it with a tumour cell population. Fibroblasts produce growth factor GF, which
		 * promotes the growth of tumour cells.
		 *
		 * After time T, we simulate treatment: a "resistant" TC population is unaffected.
		 * The "non-resistant" cell population only proliferates in the presence of GF.
		 *
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Parameters
		unsigned numberOfFibroblasts = 5;
		unsigned numberOfSensitiveTumourCells = 450;
		unsigned numberOfResistantTumourCells = 50;
		double squareDomainWidth = 20;
		const int DIM = 2;

		//double squareDomainWidthHalved = squareDomainWidth*0.5;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		// We randomly distribute nodes across the domain
		unsigned nodeNum=0;
		double x;
		double y;
		//double z;

		while(nodeNum < (numberOfFibroblasts + numberOfSensitiveTumourCells + numberOfResistantTumourCells))
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			x = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			y = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			//z = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf()-squareDomainWidthHalved);

			nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));//, z));
			nodeNum++;
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_wild_state);
		MAKE_PTR(ResistantIMOMutationState, p_resistant_state);
		MAKE_PTR(FibroblastIMOMutationState, p_fibroblast_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_fibroblast_type);

		// Create tumour cells manually
		for (unsigned i=0; i<numberOfSensitiveTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(10.0);
			p_model->SetMaxCellCycleDuration(12.0);
			p_model->SetMinCellCycleDurationByFibroblast(0.1);
			p_model->SetMaxCellCycleDurationByFibroblast(1.0);

			CellPtr p_cell(new Cell(p_wild_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					( p_model->GetMinCellCycleDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
        for (unsigned i=0; i<numberOfResistantTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(10.0);
			p_model->SetMaxCellCycleDuration(12.0);
			p_model->SetMinCellCycleDurationByFibroblast(0.1);
			p_model->SetMaxCellCycleDurationByFibroblast(1.0);

			CellPtr p_cell(new Cell(p_resistant_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					( p_model->GetMinCellCycleDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
		for (unsigned i = 0; i<numberOfFibroblasts; i++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;

            CellPtr p_cell(new Cell(p_fibroblast_state, p_model));
            p_cell->SetCellProliferativeType(p_fibroblast_type);

            cells.push_back(p_cell);
        }


		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.AddCellWriter<CellMutationStatesWriter>();
		//cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		ChastePoint<DIM> lower(0,0);
		ChastePoint<DIM> upper(squareDomainWidth,squareDomainWidth);
		//ChastePoint<3> lower(-squareDomainWidthHalved, -squareDomainWidthHalved);//, -squareDomainWidthHalved-0.1);
		//ChastePoint<3> upper(squareDomainWidthHalved, squareDomainWidthHalved);//, squareDomainWidthHalved+0.1);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		/*
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifierVariable<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("HGF");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
        */

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		//simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("IMO2017/AddingResistantTumourCellPopulation/");
		simulator.SetSamplingTimestepMultiple(12); // 12 timesteps = 6 minutes = 0.1 hours
		simulator.SetEndTime(24.0);



		// Add boundary conditions
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = unit_vector<double>(2,1);
        normal(0) = -1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc1, (&cell_population, point, normal));
        p_bc1->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        normal(0) = 0.0; normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc2, (&cell_population, point, normal));
        p_bc2->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = squareDomainWidth; point(1) = squareDomainWidth;
        normal(0) = 0.0; normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc3, (&cell_population, point, normal));
        p_bc3->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        normal(0) = 1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc4, (&cell_population, point, normal));
        p_bc4->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        //MAKE_PTR(ResistantMutationStateTrackingModifier<DIM>, p_resistanceModifier);
        //simulator.AddSimulationModifier(p_resistanceModifier);

        // Add a cell killer
        double probabiltyOfDeathInOneHour = 0.01;
        MAKE_PTR_ARGS(RandomCellKillerIgnoringFibroblastsIMO<DIM>, p_random_killer, (&cell_population, probabiltyOfDeathInOneHour));
        simulator.AddCellKiller(p_random_killer);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(15.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Label cells which are in close proximity to fibroblasts
		MAKE_PTR(MacrophageProximityLabellingModifier<DIM>, p_proximity_labeller);
		p_proximity_labeller->SetMacrophageProximityLabellingRadius(2.0);
		simulator.AddSimulationModifier(p_proximity_labeller);


		/*** To run the simulation, we call {{{Solve()}}}. ***/
		simulator.Solve();

		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

		/*** We now need to add a simulation modifier to simulate adding the treatment */



		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

    void TestAddingTreatmentRandomFibroblasts() throw(Exception)
	{
		/* We throw down fibroblasts randomly across the domain, and seed
		 * it with a tumour cell population. Fibroblasts produce growth factor GF, which
		 * promotes the growth of tumour cells.
		 *
		 * After time T, we simulate treatment: a "resistant" TC population is unaffected.
		 * The "non-resistant" cell population only proliferates in the presence of GF.
		 *
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Parameters
		unsigned numberOfFibroblasts = 250;
		unsigned numberOfSensitiveTumourCells = 250;
		unsigned numberOfResistantTumourCells = 25;
		double squareDomainWidth = 40;
		const int DIM = 2;
		double simulationDuration = 100;
		double timeToApplyTreatment = 10;
		double distanceFromFibroblastToFeelBenefit = 4.0;

		/** Treatment Parameters */
		// Pre-Treatment
		double MinCellCycleDuration_Sensitive = 0.1;
		double MaxCellCycleDuration_Sensitive = 0.2;

        double MinCellCycleDuration_Resistant = 0.2;
		double MaxCellCycleDuration_Resistant = 0.3;
		// Post-treatment
		double MinCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MaxCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MinCellCycleDurationByFibroblast_SensitivePostTreatment = 0.1;
		double MaxCellCycleDurationByFibroblast_SensitivePostTreatment = 0.2;

		std::string outputDirectory = "IMO2017/RandomDistribution/Treatment_Tis10_LargeBox/";

		//double squareDomainWidthHalved = squareDomainWidth*0.5;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		// We randomly distribute nodes across the domain
		unsigned nodeNum=0;
		double x;
		double y;
		//double z;

		while(nodeNum < (numberOfFibroblasts + numberOfSensitiveTumourCells + numberOfResistantTumourCells))
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			x = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			y = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			//z = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf()-squareDomainWidthHalved);

			nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));//, z));
			nodeNum++;
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_wild_state);
		MAKE_PTR(ResistantIMOMutationState, p_resistant_state);
		MAKE_PTR(FibroblastIMOMutationState, p_fibroblast_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_fibroblast_type);

		// Create tumour cells manually
		for (unsigned i=0; i<numberOfSensitiveTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Sensitive);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);

			CellPtr p_cell(new Cell(p_wild_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					( p_model->GetMinCellCycleDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
        for (unsigned i=0; i<numberOfResistantTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Resistant);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);

			CellPtr p_cell(new Cell(p_resistant_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_model->GetMinCellCycleDuration();
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
		for (unsigned i = 0; i<numberOfFibroblasts; i++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;

            CellPtr p_cell(new Cell(p_fibroblast_state, p_model));
            p_cell->SetCellProliferativeType(p_fibroblast_type);

            cells.push_back(p_cell);
        }


		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.AddCellWriter<CellMutationStatesWriter>();
		//cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		ChastePoint<DIM> lower(0,0);
		ChastePoint<DIM> upper(squareDomainWidth,squareDomainWidth);
		//ChastePoint<3> lower(-squareDomainWidthHalved, -squareDomainWidthHalved);//, -squareDomainWidthHalved-0.1);
		//ChastePoint<3> upper(squareDomainWidthHalved, squareDomainWidthHalved);//, squareDomainWidthHalved+0.1);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		/*
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifierVariable<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("HGF");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
        */

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		//simulator.AddSimulationModifier(p_pde_modifier);

		simulator.SetOutputDirectory(outputDirectory);
		simulator.SetSamplingTimestepMultiple(24); // 12 timesteps = 6 minutes = 0.1 hours
		simulator.SetEndTime(simulationDuration);



		// Add boundary conditions
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = unit_vector<double>(2,1);
        normal(0) = -1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc1, (&cell_population, point, normal));
        p_bc1->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        normal(0) = 0.0; normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc2, (&cell_population, point, normal));
        p_bc2->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = squareDomainWidth; point(1) = squareDomainWidth;
        normal(0) = 0.0; normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc3, (&cell_population, point, normal));
        p_bc3->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        normal(0) = 1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc4, (&cell_population, point, normal));
        p_bc4->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(AddProliferationInhibitingTreatmentModifier<DIM>, p_treatment_modifier);
        p_treatment_modifier->SetTreatmentStartTime(timeToApplyTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MinCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MaxCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MinCellCycleDurationByFibroblast_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MaxCellCycleDurationByFibroblast_SensitivePostTreatment);
        simulator.AddSimulationModifier(p_treatment_modifier);

        // Add a cell killer
        double probabiltyOfDeathInOneHour = 0.05;
        MAKE_PTR_ARGS(RandomCellKillerIgnoringFibroblastsIMO<DIM>, p_random_killer, (&cell_population, probabiltyOfDeathInOneHour));
        simulator.AddCellKiller(p_random_killer);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		//MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		//p_force->SetMeinekeSpringStiffness(15.0);
		//p_force->SetCutOffLength(1.5);
		//simulator.AddForce(p_force);

		//OR


        MAKE_PTR(RepulsionForceKeepMacrophagesStationary<DIM>, p_force);
        simulator.AddForce(p_force);

        /*
        MAKE_PTR(KeepMacrophagesStationaryLinearSpringForce<DIM>, p_force);
        p_force->SetMeinekeSpringStiffness(15.0);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);
        */




		// Label cells which are in close proximity to fibroblasts
		MAKE_PTR(MacrophageProximityLabellingModifier<DIM>, p_proximity_labeller);
		p_proximity_labeller->SetMacrophageProximityLabellingRadius(distanceFromFibroblastToFeelBenefit);
		simulator.AddSimulationModifier(p_proximity_labeller);


		/*** To run the simulation, we call {{{Solve()}}}. ***/
		simulator.Solve();

		//CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void dontTestReadFibroblastsFromData() throw(Exception)
	{
		/* We throw down fibroblasts randomly across the domain, and seed
		 * it with a tumour cell population. Fibroblasts produce growth factor GF, which
		 * promotes the growth of tumour cells.
		 *
		 * After time T, we simulate treatment: a "resistant" TC population is unaffected.
		 * The "non-resistant" cell population only proliferates in the presence of GF.
		 *
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

        ColumnDataReader reader("IMO2017/FibroblastDistributions/", "GridFibs.dat");
        std::vector<double> x_values = reader.GetValues("x");
        std::vector<double> y_values = reader.GetValues("y");

		// Parameters
		unsigned numberOfFibroblasts = 5;
		unsigned numberOfSensitiveTumourCells = 500;
		unsigned numberOfResistantTumourCells = 50;
		double squareDomainWidth = 20;
		const int DIM = 2;
		double simulationDuration = 60;
		double timeToApplyTreatment = 5;
		double distanceFromFibroblastToFeelBenefit = 3.0;

		/** Treatment Parameters */
		// Pre-Treatment
		double MinCellCycleDuration_Sensitive = 0.1;
		double MaxCellCycleDuration_Sensitive = 0.2;

        double MinCellCycleDuration_Resistant = 0.2;
		double MaxCellCycleDuration_Resistant = 0.3;
		// Post-treatment
		double MinCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MaxCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MinCellCycleDurationByFibroblast_SensitivePostTreatment = 0.1;
		double MaxCellCycleDurationByFibroblast_SensitivePostTreatment = 0.2;

		//double squareDomainWidthHalved = squareDomainWidth*0.5;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		// We randomly distribute nodes across the domain
		unsigned nodeNum=0;
		double x;
		double y;
		//double z;

		while(nodeNum < (numberOfFibroblasts + numberOfSensitiveTumourCells + numberOfResistantTumourCells))
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			x = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			y = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			//z = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf()-squareDomainWidthHalved);

			nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));//, z));
			nodeNum++;
		}

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_wild_state);
		MAKE_PTR(ResistantIMOMutationState, p_resistant_state);
		MAKE_PTR(FibroblastIMOMutationState, p_fibroblast_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_fibroblast_type);

		// Create tumour cells manually
		for (unsigned i=0; i<numberOfSensitiveTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Sensitive);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);

			CellPtr p_cell(new Cell(p_wild_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					( p_model->GetMinCellCycleDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
        for (unsigned i=0; i<numberOfResistantTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Resistant);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);

			CellPtr p_cell(new Cell(p_resistant_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_model->GetMinCellCycleDuration();
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
		for (unsigned i = 0; i<numberOfFibroblasts; i++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;

            CellPtr p_cell(new Cell(p_fibroblast_state, p_model));
            p_cell->SetCellProliferativeType(p_fibroblast_type);

            cells.push_back(p_cell);
        }


		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.AddCellWriter<CellMutationStatesWriter>();
		//cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		ChastePoint<DIM> lower(0,0);
		ChastePoint<DIM> upper(squareDomainWidth,squareDomainWidth);
		//ChastePoint<3> lower(-squareDomainWidthHalved, -squareDomainWidthHalved);//, -squareDomainWidthHalved-0.1);
		//ChastePoint<3> upper(squareDomainWidthHalved, squareDomainWidthHalved);//, squareDomainWidthHalved+0.1);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		/*
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifierVariable<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("HGF");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
        */

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		//simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("IMO2017/PreloadedFibroblasts/");
		simulator.SetSamplingTimestepMultiple(30); // 12 timesteps = 6 minutes = 0.1 hours
		simulator.SetEndTime(simulationDuration);



		// Add boundary conditions
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = unit_vector<double>(2,1);
        normal(0) = -1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc1, (&cell_population, point, normal));
        p_bc1->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        normal(0) = 0.0; normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc2, (&cell_population, point, normal));
        p_bc2->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = squareDomainWidth; point(1) = squareDomainWidth;
        normal(0) = 0.0; normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc3, (&cell_population, point, normal));
        p_bc3->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        normal(0) = 1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc4, (&cell_population, point, normal));
        p_bc4->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(AddProliferationInhibitingTreatmentModifier<DIM>, p_treatment_modifier);
        p_treatment_modifier->SetTreatmentStartTime(timeToApplyTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MinCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MaxCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MinCellCycleDurationByFibroblast_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MaxCellCycleDurationByFibroblast_SensitivePostTreatment);
        simulator.AddSimulationModifier(p_treatment_modifier);

        // Add a cell killer
        double probabiltyOfDeathInOneHour = 0.04;
        MAKE_PTR_ARGS(RandomCellKillerIgnoringFibroblastsIMO<DIM>, p_random_killer, (&cell_population, probabiltyOfDeathInOneHour));
        simulator.AddCellKiller(p_random_killer);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(15.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Label cells which are in close proximity to fibroblasts
		MAKE_PTR(MacrophageProximityLabellingModifier<DIM>, p_proximity_labeller);
		p_proximity_labeller->SetMacrophageProximityLabellingRadius(distanceFromFibroblastToFeelBenefit);
		simulator.AddSimulationModifier(p_proximity_labeller);


		/*** To run the simulation, we call {{{Solve()}}}. ***/
		simulator.Solve();

		//CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void dontTestAddingTreatmentInGrid() throw(Exception)
	{
		/* We throw down fibroblasts randomly across the domain, and seed
		 * it with a tumour cell population. Fibroblasts produce growth factor GF, which
		 * promotes the growth of tumour cells.
		 *
		 * After time T, we simulate treatment: a "resistant" TC population is unaffected.
		 * The "non-resistant" cell population only proliferates in the presence of GF.
		 *
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Parameters
		//unsigned numberOfFibroblasts = 100;
		unsigned numberOfSensitiveTumourCells = 5000;
		unsigned numberOfResistantTumourCells = 500;
		double squareDomainWidth = 40;
		const int DIM = 2;
		double simulationDuration = 2000;
		double timeToApplyTreatment = 1;
		double distanceFromFibroblastToFeelBenefit = 4.0;

		/** Treatment Parameters */
		// Pre-Treatment
		double MinCellCycleDuration_Sensitive = 0.1;
		double MaxCellCycleDuration_Sensitive = 0.2;

        double MinCellCycleDuration_Resistant = 0.2;
		double MaxCellCycleDuration_Resistant = 0.3;
		// Post-treatment
		double MinCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MaxCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MinCellCycleDurationByFibroblast_SensitivePostTreatment = 0.1;
		double MaxCellCycleDurationByFibroblast_SensitivePostTreatment = 0.2;

		//double squareDomainWidthHalved = squareDomainWidth*0.5;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		// We randomly distribute nodes across the domain
		unsigned nodeNum=0;
		double x;
		double y;
		//double z;

		while(nodeNum < (numberOfSensitiveTumourCells + numberOfResistantTumourCells))
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			x = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			y = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			//z = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf()-squareDomainWidthHalved);

			nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));//, z));
			nodeNum++;
		}

		unsigned numberOfFibroblasts = 0;
		for(int x = 5; x < squareDomainWidth; x = x + 5)
        {
            for(int y = 5; y < squareDomainWidth; y = y + 5)
            {
                nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));//, z));
                nodeNum++;
                numberOfFibroblasts++;
            }
        }

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_wild_state);
		MAKE_PTR(ResistantIMOMutationState, p_resistant_state);
		MAKE_PTR(FibroblastIMOMutationState, p_fibroblast_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_fibroblast_type);

		// Create tumour cells manually
		for (unsigned i=0; i<numberOfSensitiveTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Sensitive);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);

			CellPtr p_cell(new Cell(p_wild_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					( p_model->GetMinCellCycleDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
        for (unsigned i=0; i<numberOfResistantTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Resistant);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);

			CellPtr p_cell(new Cell(p_resistant_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_model->GetMinCellCycleDuration();
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
		for (unsigned i = 0; i<numberOfFibroblasts; i++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;

            CellPtr p_cell(new Cell(p_fibroblast_state, p_model));
            p_cell->SetCellProliferativeType(p_fibroblast_type);

            cells.push_back(p_cell);
        }


		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.AddCellWriter<CellMutationStatesWriter>();
		//cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		ChastePoint<DIM> lower(0,0);
		ChastePoint<DIM> upper(squareDomainWidth,squareDomainWidth);
		//ChastePoint<3> lower(-squareDomainWidthHalved, -squareDomainWidthHalved);//, -squareDomainWidthHalved-0.1);
		//ChastePoint<3> upper(squareDomainWidthHalved, squareDomainWidthHalved);//, squareDomainWidthHalved+0.1);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		/*
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifierVariable<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("HGF");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
        */

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		//simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("IMO2017/Grid_large/");
		simulator.SetSamplingTimestepMultiple(60); // 12 timesteps = 6 minutes = 0.1 hours
		simulator.SetEndTime(simulationDuration);



		// Add boundary conditions
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = unit_vector<double>(2,1);
        normal(0) = -1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc1, (&cell_population, point, normal));
        p_bc1->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        normal(0) = 0.0; normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc2, (&cell_population, point, normal));
        p_bc2->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = squareDomainWidth; point(1) = squareDomainWidth;
        normal(0) = 0.0; normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc3, (&cell_population, point, normal));
        p_bc3->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        normal(0) = 1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc4, (&cell_population, point, normal));
        p_bc4->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(AddProliferationInhibitingTreatmentModifier<DIM>, p_treatment_modifier);
        p_treatment_modifier->SetTreatmentStartTime(timeToApplyTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MinCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MaxCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MinCellCycleDurationByFibroblast_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MaxCellCycleDurationByFibroblast_SensitivePostTreatment);
        simulator.AddSimulationModifier(p_treatment_modifier);

        // Add a cell killer
        double probabiltyOfDeathInOneHour = 0.05;
        MAKE_PTR_ARGS(RandomCellKillerIgnoringFibroblastsIMO<DIM>, p_random_killer, (&cell_population, probabiltyOfDeathInOneHour));
        simulator.AddCellKiller(p_random_killer);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(15.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Label cells which are in close proximity to fibroblasts
		MAKE_PTR(MacrophageProximityLabellingModifier<DIM>, p_proximity_labeller);
		p_proximity_labeller->SetMacrophageProximityLabellingRadius(distanceFromFibroblastToFeelBenefit);
		simulator.AddSimulationModifier(p_proximity_labeller);


		/*** To run the simulation, we call {{{Solve()}}}. ***/
		simulator.Solve();

		//CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

    void dontTestAddingTreatmentInCircles() throw(Exception)
	{
		/* We throw down fibroblasts randomly across the domain, and seed
		 * it with a tumour cell population. Fibroblasts produce growth factor GF, which
		 * promotes the growth of tumour cells.
		 *
		 * After time T, we simulate treatment: a "resistant" TC population is unaffected.
		 * The "non-resistant" cell population only proliferates in the presence of GF.
		 *
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Parameters
		//unsigned numberOfFibroblasts = 100;
		unsigned numberOfSensitiveTumourCells = 5000;
		unsigned numberOfResistantTumourCells = 500;
		double squareDomainWidth = 40;
		const int DIM = 2;
		double simulationDuration = 2000;
		double timeToApplyTreatment = 1;
		double distanceFromFibroblastToFeelBenefit = 4.0;

		/** Treatment Parameters */
		// Pre-Treatment
		double MinCellCycleDuration_Sensitive = 0.1;
		double MaxCellCycleDuration_Sensitive = 0.2;

        double MinCellCycleDuration_Resistant = 0.2;
		double MaxCellCycleDuration_Resistant = 0.3;
		// Post-treatment
		double MinCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MaxCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MinCellCycleDurationByFibroblast_SensitivePostTreatment = 0.1;
		double MaxCellCycleDurationByFibroblast_SensitivePostTreatment = 0.2;

		//double squareDomainWidthHalved = squareDomainWidth*0.5;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		// We randomly distribute nodes across the domain
		unsigned nodeNum=0;
		double x;
		double y;
		//double z;

		while(nodeNum < (numberOfSensitiveTumourCells + numberOfResistantTumourCells))
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			x = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			y = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			//z = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf()-squareDomainWidthHalved);

			nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));//, z));
			nodeNum++;
		}

		unsigned numberOfFibroblasts = 0;
		for(int x = 10; x < squareDomainWidth; x = x + 20)
        {
            for(int y = 10; y < squareDomainWidth; y = y + 20)
            {
                for(int i = -6; i < 6; i++ )
                {
                    for(int j = -6; j < 6; j++)
                    {
                        if((i*i + j*j > 4) && (i*i + j*j < 6))
                        {
                            nodes.push_back(new Node<DIM>(nodeNum,  false,  x + i, y + j));//, z));
                            nodeNum++;
                            numberOfFibroblasts++;
                        }
                    }

                }
            }
        }

		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_wild_state);
		MAKE_PTR(ResistantIMOMutationState, p_resistant_state);
		MAKE_PTR(FibroblastIMOMutationState, p_fibroblast_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_fibroblast_type);

		// Create tumour cells manually
		for (unsigned i=0; i<numberOfSensitiveTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Sensitive);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);

			CellPtr p_cell(new Cell(p_wild_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					( p_model->GetMinCellCycleDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
        for (unsigned i=0; i<numberOfResistantTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Resistant);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);

			CellPtr p_cell(new Cell(p_resistant_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_model->GetMinCellCycleDuration();
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
		for (unsigned i = 0; i<numberOfFibroblasts; i++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;

            CellPtr p_cell(new Cell(p_fibroblast_state, p_model));
            p_cell->SetCellProliferativeType(p_fibroblast_type);

            cells.push_back(p_cell);
        }


		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.AddCellWriter<CellMutationStatesWriter>();
		//cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		ChastePoint<DIM> lower(0,0);
		ChastePoint<DIM> upper(squareDomainWidth,squareDomainWidth);
		//ChastePoint<3> lower(-squareDomainWidthHalved, -squareDomainWidthHalved);//, -squareDomainWidthHalved-0.1);
		//ChastePoint<3> upper(squareDomainWidthHalved, squareDomainWidthHalved);//, squareDomainWidthHalved+0.1);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		/*
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifierVariable<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("HGF");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
        */

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		//simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("IMO2017/Circles_large/");
		simulator.SetSamplingTimestepMultiple(12); // 12 timesteps = 6 minutes = 0.1 hours
		simulator.SetEndTime(simulationDuration);



		// Add boundary conditions
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = unit_vector<double>(2,1);
        normal(0) = -1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc1, (&cell_population, point, normal));
        p_bc1->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        normal(0) = 0.0; normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc2, (&cell_population, point, normal));
        p_bc2->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = squareDomainWidth; point(1) = squareDomainWidth;
        normal(0) = 0.0; normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc3, (&cell_population, point, normal));
        p_bc3->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        normal(0) = 1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc4, (&cell_population, point, normal));
        p_bc4->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(AddProliferationInhibitingTreatmentModifier<DIM>, p_treatment_modifier);
        p_treatment_modifier->SetTreatmentStartTime(timeToApplyTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MinCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MaxCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MinCellCycleDurationByFibroblast_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MaxCellCycleDurationByFibroblast_SensitivePostTreatment);
        simulator.AddSimulationModifier(p_treatment_modifier);

        // Add a cell killer
        double probabiltyOfDeathInOneHour = 0.05;
        MAKE_PTR_ARGS(RandomCellKillerIgnoringFibroblastsIMO<DIM>, p_random_killer, (&cell_population, probabiltyOfDeathInOneHour));
        simulator.AddCellKiller(p_random_killer);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(15.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);


		// Label cells which are in close proximity to fibroblasts
		MAKE_PTR(MacrophageProximityLabellingModifier<DIM>, p_proximity_labeller);
		p_proximity_labeller->SetMacrophageProximityLabellingRadius(distanceFromFibroblastToFeelBenefit);
		simulator.AddSimulationModifier(p_proximity_labeller);


		/*** To run the simulation, we call {{{Solve()}}}. ***/
		simulator.Solve();

		//CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}


};



#endif /* TEST_IMO17_FIBROBLASTTUMOURPATTERNS_HPP_ */

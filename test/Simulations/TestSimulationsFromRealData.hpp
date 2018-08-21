#ifndef TESTSIMULATIONSFROMREALDATA_HPP_
#define TESTSIMULATIONSFROMREALDATA_HPP_


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
#include "ExecutableSupport.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodeLocationWriter.hpp"

#include "SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

#include "ColumnDataReader.hpp"
#include "FunctionalBoundaryCondition.hpp"

#include "AveragedSourceEllipticPde.hpp"
#include "AveragedSourceEllipticPdeOxygenBelowThreshold.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticBoxDomainPdeModifierVariableTimestep.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "ApoptoticCellKiller.hpp"

#include "DiffusionForceChooseD.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "NodeLocationWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "BoundaryCellWriter.hpp"
#include "BoundaryNodeWriter.hpp"

#include "ExternalPressureForceOnConcaveHull.hpp"
#include "AddMacrophagesAtSpecifiedTimeModifier.hpp"

#include "CuboidPeriodicBoundaryCondition.hpp"

#include "OutputOnlyMacrophageSummaryStatisticsModifierWithCSF1.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages.hpp"
#include "UniformCellCycleModelWithQuiescence.hpp"
#include "ChemotacticForceCSF1.hpp"
#include "OutputOnlyMacrophageSummaryStatisticsModifierWithCSF1.hpp"
#include "EllipticBoxDomainPdeModifierSetOxygenFromFile.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>


double bc_func(const ChastePoint<2>& p)
{

	int x = std::floor(p[0]);
	int y = std::floor(p[1]);
	std::cout << "Applying BCs to point at x = " << std::to_string(x) << ", y = " << std::to_string(y) << std::endl;


	ifstream inFile;

	inFile.open("/mi/share/scratch/bull/ChasteStuff/JoshuaBull/Data/22680-08_ROI_33_predictedO2Mask.csv");
	if (!inFile)
	{
		std::cout << "Unable to open file A";
		exit(1); // terminate with error
	}

	std::string value;
	std::string s;
	while (inFile)
	{
		// Keep getline-ing until we reach row y
		for(int i = 0; i < y; i++)
		{
			std::getline(inFile, value,'\n');
		}
		std::getline(inFile, value,'\n');

		std::istringstream ss( value );
		// Now find element x of that row
		for(int i = 0; i < x; i++)
		{
			std::getline(ss, s,',');
		}
		std::getline(ss, s,',');
		std::cout << s << std::endl;
		break;
	}

	inFile.close();

	return std::stod(s);
}




class TestSimulationsFromRealData : public AbstractCellBasedTestSuite
{
public:


	void TestImportHistologyDataForTumourEnvironment() throw(Exception)
	{
		const int DIM = 2;
		double squareMeshWidth = 100;

		double hypoxicConcentration = 0.3;
		double quiescentConcentration = 0.5;
		double criticalHypoxicDuration = 8;
		double averageCellCycleLength = 16;
		double apoptosisDuration = 48;


		// Prepare matrix of stroma values: 1 = stroma, 0 = tumour (PanCK)
		ifstream inFile;
		std::cout << "Reading data" << std::endl;

		inFile.open("/mi/share/scratch/bull/ChasteStuff/JoshuaBull/Data/22680-08_ROI_33_stromaLocations.csv");
		if (!inFile)
		{
			std::cout << "Unable to open file B";
			exit(1); // terminate with error
		}

		std::string value;
		std::string s;

		std::vector<std::vector<double> > stromaValues;
		while (inFile)
		{
			for(int i = 0; i < 100; i++)
			{
				std::vector<double> newRow;
				std::getline(inFile, value,'\n');
				std::istringstream ss( value );

				for(int j = 0; j < 100; j++)
				{
					std::getline(ss, s,',');
					newRow.push_back(std::stod(s));
					//std::cout << std::stod(s) << std::endl;
				}
				stromaValues.push_back(newRow);
			}
			break;
		}

		inFile.close();

		std::cout << "Data is read successfully" << std::endl;
		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;

		double simDuration = 300;

		// Add dense grid of nodes across domain
		unsigned nodeNum=0;

		for (double x=0; x<squareMeshWidth; x=x+0.8)
		{
			for (double y=0; y<squareMeshWidth; y=y+0.8)
			{
				nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));
				nodeNum++;
			}
		}

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		// All cells follow same cell cycle model, but stromal cells differentiate slower
		for (unsigned i=0; i<nodeNum; i++)
		{
			std::cout << "nodeNum = " << std::to_string(i) << " of "<< std::to_string(nodeNum) << std::endl;
			c_vector<double,DIM> location = nodes[i]->rGetLocation();
			int x = std::floor(location[0]);
			int y = std::floor(location[1]);

			bool isStromal = false;
			if(stromaValues[y][x] > 0)
			{

				isStromal = true;
			}
			std::cout << isStromal << std::endl;

			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*0.75);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*1.25);
			if(isStromal)
			{
				p_model->SetMinCellCycleDuration(32*0.75);
				p_model->SetMaxCellCycleDuration(32*1.25);
			}
			p_model->SetHypoxicConcentration(hypoxicConcentration);
			p_model->SetQuiescentConcentration(quiescentConcentration);
			p_model->SetCriticalHypoxicDuration(criticalHypoxicDuration);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->GetCellData()->SetItem("csf1", 0);
			p_cell->SetApoptosisTime(apoptosisDuration); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.75);
			if(isStromal)
			{
				// Stromal cells have 32 hour cell cycle
				birthTime = - RandomNumberGenerator::Instance()->ranf() * (32*0.75);
			}
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
		// Make PDE (Oxygen)
		double diffusionCoefficient = 1;
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_pde, (cell_population, 0, diffusionCoefficient));
		MAKE_PTR_ARGS(FunctionalBoundaryCondition<DIM>, p_functional_bc, (&bc_func));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Create a ChasteCuboid on which to base the oxygen finite element mesh
		ChastePoint<DIM> lower(0, 0);
		ChastePoint<DIM> upper(squareMeshWidth, squareMeshWidth);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		// We massively cheat, by specifying an elliptic PDE but specifying the O2 at every point of the mesh as a BC. Then, we never update the grid...
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		int updateIntervalForPdeInTimesteps = 120;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierSetOxygenFromFile<DIM>, p_pde_modifier, (p_pde, p_functional_bc, is_neumann_bc, p_cuboid, 1.0));
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(false); //was false

		//		// Make PDE (CSF1)
		//		// Assume CSF1 diffuses slower than o2 (0.2)
		//		MAKE_PTR_ARGS(AveragedSourceEllipticPdeOxygenBelowThreshold<DIM>, p_pdeCSF, (cell_population, 0.01,0.2,hypoxicConcentration));
		//		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bcCSF, (0.0));
		//
		//		bool is_neumann_bcCSF = false; // Dirichlet BCs
		//		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		//		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboidCSF, (lower, upper));
		//		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboidCSF, 1.0));
		//		p_pde_modifierCSF->SetTimestepInterval(updateIntervalForPdeInTimesteps);
		//		p_pde_modifierCSF->SetDependentVariableName("csf1");
		//		p_pde_modifierCSF->SetOutputGradient(true);
		//		p_pde_modifierCSF->SetBcsOnBoxBoundary(true); //was false


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		//		simulator.AddSimulationModifier(p_pde_modifierCSF);
		std::stringstream output_directory;
		output_directory << "SimulationsFromRealData/";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(60); // Every 4 hours
		simulator.SetEndTime(simDuration);

		//		// Add macrophages at set time
		//		MAKE_PTR(AddMacrophagesAtSpecifiedTimeModifier<DIM>, p_addMacs_modifier);
		//		p_addMacs_modifier->SetNumberOfMacrophagesToAdd(100);
		//		p_addMacs_modifier->SetTimeToAddMacrophages(timeToAddMacrophages);
		//		//p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(250);
		//		simulator.AddSimulationModifier(p_addMacs_modifier);


		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifierWithCSF1<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescentConcentration);
		p_macStats_modifier->SetHypoxiaLevel(hypoxicConcentration);
		p_macStats_modifier->SetOutputFrequencyInHours(0.5); // Every 30 minutes
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.AddSimulationModifier(p_macStats_modifier);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);

		//		// Add Chemotaxis
		//		MAKE_PTR(ChemotacticForceCSF1<DIM>, p_chemotactic_force);
		//		p_chemotactic_force->SetChemotaxisSensitivity(chemotaxisCoefficient);
		//		simulator.AddForce(p_chemotactic_force);


		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(45.0);
		p_force->SetCutOffLength(1.5);
		p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
		p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
		simulator.AddForce(p_force);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);



		// Add boundary conditions - keep things in the box
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

		point(0) = squareMeshWidth; point(1) = squareMeshWidth;
		normal(0) = 0.0; normal(1) = 1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc3, (&cell_population, point, normal));
		p_bc3->SetUseJiggledNodesOnPlane(true);
		simulator.AddCellPopulationBoundaryCondition(p_bc3);
		normal(0) = 1.0; normal(1) = 0.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc4, (&cell_population, point, normal));
		p_bc4->SetUseJiggledNodesOnPlane(true);
		simulator.AddCellPopulationBoundaryCondition(p_bc4);










		//		MAKE_PTR(ExternalPressureForceOnConcaveHull<DIM>, p_pressure);
		//		p_pressure->SetPressure(5);
		//		simulator.AddForce(p_pressure);

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


#endif /* TESTSIMULATIONSFROMREALDATA_HPP_ */

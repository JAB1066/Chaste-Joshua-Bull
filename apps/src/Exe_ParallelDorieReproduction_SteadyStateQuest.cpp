/*

Copyright (c) 2005-2016, University of Oxford.
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

#include <cxxtest/TestSuite.h>
//#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
//#include "PetscSetupAndFinalize.hpp"

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
#include "UniformCellCycleModelWithQuiescence.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

#include "AveragedSourceEllipticPde.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "EllipticBoxDomainPdeModifierVariableTimestep.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "ApoptoticCellKiller.hpp"

#include "DiffusionForceChooseD.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"

#include "CuboidPeriodicBoundaryCondition.hpp"

#include "OutputOnlyMacrophageSummaryStatisticsModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



/*
 * Prototype functions
 */
void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(std::string id_string, double averageCellCycleLength, double hypoxicConcentration, double quiescentConcentration, double criticalHypoxicDuration);
void OutputToConsole(std::string id_string, std::string leading);

int main(int argc, char *argv[])
{
	// This sets up PETSc and prints out copyright information, etc.
	//ExecutableSupport::StandardStartup(&argc, &argv);
	ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

	// Define command line options
	boost::program_options::options_description general_options("This is a Chaste executable.\n");
	general_options.add_options()
                    										("help", "produce help message")
															("ID", boost::program_options::value<std::string>(),"ID string for the simulation")
															("ACCL", boost::program_options::value<double>()->default_value(13.0),"averageCellCycleLength: cell cycle takes 1 hour either side of this duration")
															("HC", boost::program_options::value<double>()->default_value(0.15),"hypoxicConcentration: O2 conc at which cell becomes hypoxic")
															("QC", boost::program_options::value<double>()->default_value(0.35),"quiescentConcentration: O2 conc at which cell becomes quiescent")
															("CHD", boost::program_options::value<double>()->default_value(8.0),"criticalHypoxicDuration: duration cell must be hypoxic for before death");

	// Define parse command line into variables_map
	boost::program_options::variables_map variables_map;
	boost::program_options::store(parse_command_line(argc, argv, general_options), variables_map);

	// Print help message if wanted
	if (variables_map.count("help"))
	{
		std::cout << setprecision(3) << general_options << "\n";
		std::cout << general_options << "\n";
		return 1;
	}

	// Get ID and name from command line
	std::string id_string = variables_map["ID"].as<std::string>();
	double averageCellCycleLength = variables_map["ACCL"].as<double>();
	double hypoxicConcentration = variables_map["HC"].as<double>();
	double quiescentConcentration = variables_map["QC"].as<double>();
	double criticalHypoxicDuration = variables_map["CHD"].as<double>();

	OutputToConsole(id_string, "Started");

	SetupSingletons();
	SetupAndRunSimulation(id_string,averageCellCycleLength,hypoxicConcentration,quiescentConcentration,criticalHypoxicDuration);
	DestroySingletons();

	OutputToConsole(id_string, "Completed");
}

void SetupSingletons()
{
	// Set up what the test suite would do
	SimulationTime::Instance()->SetStartTime(0.0);

	// Reseed with 0 for same random numbers each time, or time(NULL) or simulation_id to change each realisation
	RandomNumberGenerator::Instance()->Reseed(time(NULL));
	CellPropertyRegistry::Instance()->Clear();
	CellId::ResetMaxCellId();
}

void DestroySingletons()
{
	// This is from the tearDown method of the test suite
	SimulationTime::Destroy();
	RandomNumberGenerator::Destroy();
	CellPropertyRegistry::Instance()->Clear();
}

void OutputToConsole(std::string id_string, std::string leading)
{
	// Compose the message
	std::stringstream message;
	message << leading << " simulation with ID string " << id_string << std::endl;

	// Send it to the console
	std::cout << message.str() << std::flush;
}

void SetupAndRunSimulation(std::string id_string, double averageCellCycleLength, double hypoxicConcentration, double quiescentConcentration, double criticalHypoxicDuration)
{
	{
		/*
		 * In this test, "macrophages" have no macrophage properties - they are completely inert,
		 * with the exception of when a degree of Brownian Motion applies to the entire simulation.
		 *
		 * As a result, we can interpret them as inert "microspheres" a la Dorie 1982
		 *
		 */

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
		double startTime = 240.0;
		double simDuration = 300.0;

		OffLatticeSimulation<DIM>* p_archivedSimulator = CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Load("BenchmarkTumour/2D_120hours", startTime);

		NodeBasedCellPopulation<DIM>& r_archivedPopulation = dynamic_cast<NodeBasedCellPopulation<DIM> &>(p_archivedSimulator->rGetCellPopulation());

		std::list<CellPtr> archived_cells = r_archivedPopulation.rGetCells();
		std::list<CellPtr>::iterator it;
		int nodeNum = 0;

		for (it = archived_cells.begin(); it != archived_cells.end(); it++)
		{
			// Give each cell a pointer to the property registry (we have taken ownership in this constructor)
			c_vector<double, DIM> location = r_archivedPopulation.GetLocationOfCellCentre(*it);
			nodes.push_back(new Node<DIM>(nodeNum,  false,  location[0], location[1], location[2]));
			double o2conc = (*it)->GetCellData()->GetItem("oxygen");
			double birthTime = (*it)->GetBirthTime();


			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength-1.0);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength+1.0);
			p_model->SetHypoxicConcentration(hypoxicConcentration);
			p_model->SetQuiescentConcentration(quiescentConcentration);
			p_model->SetCriticalHypoxicDuration(criticalHypoxicDuration);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", o2conc);
			p_cell->SetApoptosisTime(1.0); // Apoptosis time in hours - how long before a cell is removed?
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

			nodeNum++;
		}


//		/*
//		 * Now we add our macrophages
//		 */
//		unsigned firstMacrophageNode = nodes.size();
//
//		// Macrophage Nodes - in a shell at fixed radius
//		// Number of Macrophages remains constant throughout simulation
//		unsigned numMacrophages=100;
//
//		std::vector<double> xVec;
//		std::vector<double> yVec;
//		std::vector<double> zVec;
//		unsigned macCount=0;
//		double random_x;
//		double random_y;
//		double zCoord;
//
//		// Add macrophages at random points at edge of spheroid
//		double macrophageSphereRadius = 9.5;
//		while(macCount < numMacrophages)
//		{
//			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
//			random_x = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
//			random_y = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
//
//			// Choose sign of z coord randomly
//			double randNum = RandomNumberGenerator::Instance()->ranf();
//			if(randNum < 0.5)
//			{
//				zCoord = sqrt(abs(pow(macrophageSphereRadius,2) - pow(random_x,2) - pow(random_y,2)));
//			}
//			else
//			{
//				zCoord = -sqrt(abs(pow(macrophageSphereRadius,2) - pow(random_x,2) - pow(random_y,2)));
//			}
//			macCount++;
//			nodes.push_back(new Node<3>(nodeNum,  false,  random_x, random_y, zCoord));
//			nodeNum++;
//		}


		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


//		// Now create macrophages
//		for (unsigned i=firstMacrophageNode; i<nodes.size(); i++)
//		{
//			// Make Macrophage
//			NoCellCycleModel* p_model2 = new NoCellCycleModel;
//			p_model2->SetDimension(3);
//
//			CellPtr p_cell(new Cell(p_state, p_model2));
//			p_cell->SetCellProliferativeType(p_macrophage_type);
//
//			cells.push_back(p_cell);
//		}

		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		//cell_population.AddPopulationWriter<NodeLocationWriter>();
		// Write summary statistic files
		//cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
		//cell_population.AddCellWriter<CellProliferativeTypesWriter>();

		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
		//if(DIM == 3)
		//{
		//ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		//ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		//}
		//if(DIM == 2)
		//{
		ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		//}
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
		int updateIntervalForPdeInTimesteps = 60;
		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);

		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<DIM>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescentConcentration);
		p_macStats_modifier->SetHypoxiaLevel(hypoxicConcentration);
		p_macStats_modifier->SetOutputFrequencyInTimesteps(60); // Every 30 minutes


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.AddSimulationModifier(p_macStats_modifier);
		std::stringstream output_directory;
		output_directory << "ParameterSweeps/DorieReproduction_SteadyStateQuest_2D/" << id_string << "/";
		simulator.SetOutputDirectory(output_directory.str());
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(UINT_MAX); // One visualisation every 30 minutes...
		simulator.SetEndTime(startTime+simDuration); // ...for 96 hours

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<DIM>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

//		// Add Brownian motion for all cells
//		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
//		p_diffusion_force->SetDiffusionScalingConstant(diffusionCoefficient);
//		simulator.AddForce(p_diffusion_force);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
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

}

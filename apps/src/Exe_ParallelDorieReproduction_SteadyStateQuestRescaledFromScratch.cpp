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
#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



/*
 * Prototype functions
 */
void SetupSingletons(int seed);
void DestroySingletons();
void SetupAndRunSimulation(std::string id_string, double iterationNumber, double averageCellCycleLength, double hypoxicConcentration, double quiescentConcentration, double criticalHypoxicDuration);
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
															("IT", boost::program_options::value<double>(),"Iteration Number")
															("ACCL", boost::program_options::value<double>()->default_value(20.0),"averageCellCycleLength: cell cycle takes 1 hour either side of this duration")
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
	double iterationNumber = variables_map["IT"].as<double>();
	double averageCellCycleLength = variables_map["ACCL"].as<double>();
	double hypoxicConcentration = variables_map["HC"].as<double>();
	double quiescentConcentration = variables_map["QC"].as<double>();
	double criticalHypoxicDuration = variables_map["CHD"].as<double>();


	OutputToConsole(id_string, "Started");

	int seed = std::stoi(id_string);
	SetupSingletons(seed);
	SetupAndRunSimulation(id_string,iterationNumber,averageCellCycleLength,hypoxicConcentration,quiescentConcentration,criticalHypoxicDuration);
	DestroySingletons();

	OutputToConsole(id_string, "Completed");
}

void SetupSingletons(int seed)
{
	// Set up what the test suite would do
	SimulationTime::Instance()->SetStartTime(0.0);

	RandomNumberGenerator::Instance()->Reseed(seed);
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

void SetupAndRunSimulation(std::string id_string, double iterationNumber, double averageCellCycleLength, double hypoxicConcentration, double quiescentConcentration, double criticalHypoxicDuration)
{
	{
		/*
		 * In this test, "macrophages" have no macrophage properties - they are completely inert,
		 * with the exception of when a degree of Brownian Motion applies to the entire simulation.
		 *
		 * As a result, we can interpret them as inert "microspheres" a la Dorie 1982
		 *
		 */

		/* We rescale our time units to effectively increase the timestep
		 */


		const int DIM = 3;
		int timestepsPerHour = 120; // Reinterpret length of a timestep
		int visualisationOutputFrequencyPerHour = 0;

		double chasteDefaultTimestepInHours = 1.0/120.0;
		double intendedNewTimestepInHours = 1.0/timestepsPerHour;
		double timeRescalingConstant = chasteDefaultTimestepInHours/intendedNewTimestepInHours;
		double simDuration = 200.0*timeRescalingConstant;

		if(hypoxicConcentration > quiescentConcentration)
		{
			std::cout << "Simulation ID " << id_string << " terminated predicatbly" << std::endl << std::flush;
			simDuration = 0;
		}



		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 40;
		
		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 3;
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

		double startTime = 0;


		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

		// Create tumour cells manually
		for (unsigned i=0; i<nodeNum; i++)
		{
			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*timeRescalingConstant*0.9);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*timeRescalingConstant*1.1);
			p_model->SetHypoxicConcentration(hypoxicConcentration);
			p_model->SetQuiescentConcentration(quiescentConcentration);
			p_model->SetCriticalHypoxicDuration(criticalHypoxicDuration*timeRescalingConstant);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->SetApoptosisTime(48.0*timeRescalingConstant); // Apoptosis time in hours - how long before a cell is removed?

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
								(averageCellCycleLength*timeRescalingConstant*0.9);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}

		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(6.0);//Set big movement threshold

		// Make PDE (Oxygen)
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_pde, (cell_population, -0.03));//was -0.03
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_bc, (1.0));
		bool is_neumann_bc = false; // Dirichlet BCs

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
		//if(DIM == 3)
		//{
		ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		//}
		//if(DIM == 2)
		//{
		//ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		//ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
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
		p_macStats_modifier->SetOutputFrequencyInTimesteps(60*timeRescalingConstant); // Every 30 minutes


		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.AddSimulationModifier(p_macStats_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/25Jan_DorieParamSweep/" << id_string << "/";
		simulator.SetOutputDirectory(output_directory.str());
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		if(visualisationOutputFrequencyPerHour > 0)
		{
			simulator.SetSamplingTimestepMultiple(timestepsPerHour/visualisationOutputFrequencyPerHour); // One visualisation every 30 minutes...
		}
		else if(visualisationOutputFrequencyPerHour == 0)
		{
			simulator.SetSamplingTimestepMultiple(UINT_MAX);
		}
		simulator.SetEndTime(startTime+simDuration); // ...for 96 hours

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<DIM>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);

//		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
//		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
//		p_force->SetMeinekeSpringStiffness(5.0);
//		p_force->SetCutOffLength(1.5);
//		simulator.AddForce(p_force);
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

		//CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

}

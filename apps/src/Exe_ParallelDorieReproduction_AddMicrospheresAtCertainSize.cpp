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
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
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
#include "NodeVelocityWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellRadiusWriter.hpp"

#include "CuboidPeriodicBoundaryCondition.hpp"

#include "OutputOnlyMacrophageSummaryStatisticsModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis.hpp"
#include "UniformCellCycleModelWithQuiescence.hpp"

#include "AddMacrophagesAtSpecifiedSpheroidSizeModifier.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



/*
 * Prototype functions
 */
void SetupSingletons(unsigned seed);
void DestroySingletons();
void SetupAndRunSimulation(std::string id_string, double iterationNumber, double averageCellCycleLength, double hypoxicConcentration, double quiescentConcentration, double criticalHypoxicDuration, double oxygenConsumptionRate);
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
																			("CHD", boost::program_options::value<double>()->default_value(8.0),"criticalHypoxicDuration: duration cell must be hypoxic for before death")
																			("OCR", boost::program_options::value<double>()->default_value(0.03),"oxygenConsumptionRate: amount of O2 consumed by each tumour cell");

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
	double oxygenConsumptionRate = variables_map["OCR"].as<double>();

	// Hack:
	if (quiescentConcentration == 0.4)
	{
		hypoxicConcentration = 0.1;
	}

	OutputToConsole(id_string, "Started");

	unsigned seed = static_cast<unsigned>(std::stoi(id_string));

	SetupSingletons(seed);
	SetupAndRunSimulation(id_string,iterationNumber,averageCellCycleLength,hypoxicConcentration,quiescentConcentration,criticalHypoxicDuration,oxygenConsumptionRate);
	DestroySingletons();

	OutputToConsole(id_string, "Completed");
}

void SetupSingletons(unsigned seed)
{
	// Set up what the test suite would do
	SimulationTime::Instance()->SetStartTime(0.0);
	//RandomNumberGenerator::Instance()->Reseed(time(NULL));
	//std::stringstream message;
	//message << "Reseeding with seed " << std::to_string(seed) << std::endl;
	//std::cout << message.str() << std::flush;
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

void SetupAndRunSimulation(std::string id_string, double iterationNumber, double averageCellCycleLength, double hypoxicConcentration, double quiescentConcentration, double criticalHypoxicDuration, double oxygenConsumptionRate)
{
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

		double simDuration = 1000.0;


		if(hypoxicConcentration > quiescentConcentration)
		{
			std::cout << "Simulation ID " << id_string << " terminated predictably" << std::endl << std::flush;
			simDuration = 0;
		}

		// Add some nodes

//		double cubeDomainDistanceToBoundary = 40;

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		// I'm suspicious of the RNG, so we generate a different number of random numbers at the start of each simulation to ensure that it starts in a different place.
		//std::srand(iterationNumber);
		unsigned seed = static_cast<unsigned>(std::stoi(id_string));
		RandomNumberGenerator::Instance()->Reseed(seed);



		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 7.0;
		for (double x=-initialRadius; x<initialRadius+1; x++)
		{
			for (double y=-initialRadius; y<initialRadius+1; y++)
			{
				for (double z=-initialRadius; z<initialRadius+1; z++)
				{
					if(pow(x,2) + pow(y,2) + pow(z,2) < pow(initialRadius,2))
					{
						nodes.push_back(new Node<3>(nodeNum,  false,  x, y, z));
						nodeNum++;
					}
				}
			}
		}

		for (unsigned i=0; i<nodeNum; i++)
		{
			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
					(averageCellCycleLength*0.9);

			UniformCellCycleModelWithQuiescence* p_model = new UniformCellCycleModelWithQuiescence;
			p_model->SetDimension(3);
			p_model->SetMinCellCycleDuration(averageCellCycleLength*0.9);
			p_model->SetMaxCellCycleDuration(averageCellCycleLength*1.1);
			p_model->SetHypoxicConcentration(hypoxicConcentration);
			p_model->SetQuiescentConcentration(quiescentConcentration);
			p_model->SetCriticalHypoxicDuration(criticalHypoxicDuration);
			//p_model->SetBirthTime(birthTime);

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1.0);
			p_cell->SetApoptosisTime(24); // Apoptosis time in hours - how long before a cell is removed?
			p_cell->SetCellCycleModel(p_model);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);
		}



		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(1000.0);//Set big movement threshold
		cell_population.AddPopulationWriter<NodeVelocityWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		//cell_population.AddCellWriter<CellRadiusWriter>();

//		// Make PDE (Oxygen)
//		MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -oxygenConsumptionRate));//was -0.03
//		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
//		bool is_neumann_bc = false; // Dirichlet BCs

//		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
//		c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
//		ChastePoint<3> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
//		ChastePoint<3> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
//		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));
//
//		// Create a PDE modifier and set the name of the dependent variable in the PDE
//		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
//		p_pde_modifier->SetDependentVariableName("oxygen");
//		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
//		int updateIntervalForPdeInTimesteps = 60;
//		p_pde_modifier->SetTimestepInterval(updateIntervalForPdeInTimesteps);
//
		// Create an output modifier - set SamplingTimestepMultiple below to UINT_MAX
		MAKE_PTR(OutputOnlyMacrophageSummaryStatisticsModifier<3>, p_macStats_modifier);
		p_macStats_modifier->SetQuiescenceLevel(quiescentConcentration);
		p_macStats_modifier->SetHypoxiaLevel(hypoxicConcentration);
		p_macStats_modifier->SetOutputFrequencyInHours(0.5); // Every 30 minutes

		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(CellwiseSourceEllipticPde<3>, p_pde, (cell_population, -oxygenConsumptionRate));
		MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
		bool is_neumann_bc = false;

		MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc));
		p_pde_modifier->SetDependentVariableName("oxygen");

//		// Add macrophages at set time
//		MAKE_PTR(AddMacrophagesAtSpecifiedSpheroidSizeModifier<3>, p_addMacs_modifier);
//		p_addMacs_modifier->SetNumberOfMacrophagesToAdd(100);
//		p_addMacs_modifier->SetNumberOfTumourCellsBeforeMacrophagesAreAdded(11000);
//		p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(100);

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
//		simulator.AddSimulationModifier(p_addMacs_modifier);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.AddSimulationModifier(p_macStats_modifier);
		std::stringstream output_directory;
		output_directory << "Dorie1982/ParamSweeps/3Mar_SteadyState/" << id_string << "/";
		p_macStats_modifier->SetOutputDirectory(output_directory.str());
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(120); // One visualisation every 30 minutes...
		simulator.SetEndTime(simDuration); // ...for 96 hours
//
//		// Add periodic boundary conditions for cells
//		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_periodic_boundary_condition,(&cell_population,lower,upper));
//		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.01);
		simulator.AddForce(p_diffusion_force);

		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosis<3>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		p_force->SetHeterotypicLabelledSpringConstantMultiplier(0.5);
		p_force->SetHomotypicLabelledSpringConstantMultiplier(0.3);
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

}

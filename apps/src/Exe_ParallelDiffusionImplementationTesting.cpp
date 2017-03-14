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

#include "CellId.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "ExecutableSupport.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines a node-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"

// Program option includes for handling command line arguments
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>

#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

#include "ApoptoticCellKiller.hpp"

#include "DiffusionForceChooseD.hpp"
#include "DiffusionForceFixedSigmaRandomUnitDirection.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "NodeLocationWriter.hpp"

#include <fstream>
#include <iostream>


/*
 * Prototype functions
 */
void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(std::string id_string, double diffusionCoefficient, double simulationDuration, double iterationNumber, double dampingConstant);
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
											("DC", boost::program_options::value<double>()->default_value(0.0),"diffusionCoefficient: amount of random motion by cell")
											("SD", boost::program_options::value<double>()->default_value(0.0),"simulationDuration: time to run (hours)")
											("IN", boost::program_options::value<double>()->default_value(0),"Iteration Number: Iterator not used in running of simulation")
											("Damp", boost::program_options::value<double>()->default_value(1),"Damping Constant: Constant for damping");




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
	double diffusionCoefficient = variables_map["DC"].as<double>();
	double simulationDuration = variables_map["SD"].as<double>();
	double iterationNumber = variables_map["IN"].as<double>();
	double dampingConstant = variables_map["Damp"].as<double>();

	OutputToConsole(id_string, "Started");

	SetupSingletons();
	SetupAndRunSimulation(id_string,diffusionCoefficient,simulationDuration,iterationNumber,dampingConstant);
	DestroySingletons();

	OutputToConsole(id_string, "Completed");
}

void SetupSingletons()
{
	// Set up what the test suite would do
	SimulationTime::Instance()->SetStartTime(0.0);

	// Reseed with 0 for same random numbers each time, or time(NULL) or simulation_id to change each realisation
	int seed = time(NULL);
	RandomNumberGenerator::Instance()->Reseed(seed);//1485955195
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


void SetupAndRunSimulation(std::string id_string, double diffusionCoefficient, double simulationDuration, double iterationNumber, double dampingConstant)
{

	// Test macrophage random motion and log location after time t
	/* We vary parameters:
	 *
	 * diffusionCoefficient: coefficient for random macrophage motion law (Fickian)
	 * simulationDuration: length of simulation in hours
	 * iterationNumber: e.g, simulation 4 of a given n for each parameter set (specified in python) - not used in simulation, just as an iterator
	 * dampingConstant: constant nu on LHS of nu*dx/dt = Forces
	 *
	 */
	const int dimensions = 3;

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
	//p_cell->GetCellData()->SetItem("csf1", 0.0);

	cells.push_back(p_cell);

	// Make cell population (2D)
	NodeBasedCellPopulation<dimensions> cell_population(mesh, cells);
	//cell_population.AddPopulationWriter<NodeLocationWriter>();
	//cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files

	cell_population.SetDampingConstantNormal(dampingConstant);
	cell_population.SetDampingConstantMutant(dampingConstant);


	/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
	 * (this time with dimension 2) and set the output directory, output multiple and end time. */
	OffLatticeSimulation<dimensions> simulator(cell_population);

	// Create and set an output directory that is different for each simulation
	std::stringstream output_directory;
	output_directory << "AddingMacrophages/DiffusionImplementationTesting/SingleMacrophage_DistanceInSetTime/" << id_string;
	simulator.SetOutputDirectory(output_directory.str());
	simulator.SetSamplingTimestepMultiple(12);
	simulator.SetEndTime(simulationDuration);



	//MAKE_PTR(DiffusionForceFixedSigmaRandomUnitDirection<dimensions>, p_diffusion_force);
	//p_diffusion_force->SetScalingConstant(diffusionCoefficient);
	//simulator.AddForce(p_diffusion_force);

	MAKE_PTR(DiffusionForceChooseD<dimensions>, p_diffusion_force);
	p_diffusion_force->SetDiffusionScalingConstant(diffusionCoefficient);
	simulator.AddForce(p_diffusion_force);


	/* To run the simulation, we call {{{Solve()}}}. */
	try
	{
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
	}
	catch(Exception& e)
	{
		throw(e);//std::cout << error << "\n";
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
						<< "1" << ",";
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
	}
	catch(...)
		{
			std::cout << "Non-exception Error caught for ID " << id_string << "\n";
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
							<< "1" << ",";
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
		}


	/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
	for (unsigned i=0; i<nodes.size(); i++)
	{
		delete nodes[i];
	}

}

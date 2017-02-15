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
#include "UniformCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "ExecutableSupport.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines a node-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"
// Write VTU files
#include "VoronoiDataWriter.hpp"
#include "RandomCellKiller.hpp"

// Program option includes for handling command line arguments
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>

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

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages.hpp"
#include "MacrophageProximityLabellingModifier.hpp"
#include "CheckIfMacrophageAtEndModifier.hpp"

#include "NodeLocationWriter.hpp"


/*
 * Prototype functions
 */
void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(std::string id_string, double lengthOfDomain, double chemotaxisSensitivity, double temperatureKelvin, double iterationNumber);
void OutputToConsole(std::string id_string, std::string leading);
double bc_func(const ChastePoint<2>& p);
double bc_func2(const ChastePoint<2>& p);

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
											("LD", boost::program_options::value<double>()->default_value(10.0),"lengthOfDomain: extent of domain in y direction")
											("CS", boost::program_options::value<double>()->default_value(1.0),"chemotaxisSensitivity: Sensitivity of MP to chemotactic gradient")
											("TK", boost::program_options::value<double>()->default_value(0.0),"Temperature in Kelvin: For the diffusion rate")
											("IN", boost::program_options::value<double>()->default_value(0),"Iteration Number: Iterator not used in running of simulation");




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
	double lengthOfDomain = variables_map["LD"].as<double>();
	double chemotaxisSensitivity = variables_map["CS"].as<double>();
	double temperatureKelvin = variables_map["TK"].as<double>();
	double iterationNumber = variables_map["IN"].as<double>();

	OutputToConsole(id_string, "Started");

	SetupSingletons();
	SetupAndRunSimulation(id_string,lengthOfDomain,chemotaxisSensitivity,temperatureKelvin,iterationNumber);
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

double bc_func(const ChastePoint<2>& p)
{
	return p[1]; // Return y value on boundary
}

double bc_func2(const ChastePoint<2>& p)
{
	return p[1]*p[1]; // Return y^2 value on boundary
}

void SetupAndRunSimulation(std::string id_string, double lengthOfDomain, double chemotaxisSensitivity, double temperatureKelvin, double iterationNumber)
{

	// Test macrophage on constant chemotactic gradient and log time to traverse.
	/* We vary parameters:
	 *
	 * lengthOfDomain: length of domain which macrophage must traverse
	 * chemotaxisSensitivity: Sensitivity of MP to chemotactic gradient
	 * temperatureKelvin: Temperature of random motion term
	 * iterationNumber: e.g, simulation 4 of a given n for each parameter set (specified in python) - not used in simulation, just as an iterator
	 *
	 */


	// Generate Mesh:
	// Make Vector
	std::vector<Node<2>*> nodes;

	// Add macrophage node at origin
	unsigned nodeNum=0;

	// Add macrophage node
	nodes.push_back(new Node<2>(nodeNum, false, 0, 0.01));

	NodesOnlyMesh<2> mesh;
	// Cut off length: 1.5 cell radii
	mesh.ConstructNodesWithoutMesh(nodes, 1.5);

	// Make cell pointers
	std::vector<CellPtr> cells;

	MAKE_PTR(WildTypeCellMutationState, p_state);
	MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

	// Make Macrophage
	NoCellCycleModel* p_model = new NoCellCycleModel;
	p_model->SetDimension(2);

	CellPtr p_cell(new Cell(p_state, p_model));
	p_cell->SetCellProliferativeType(p_macrophage_type);
	p_cell->GetCellData()->SetItem("csf1", 0.0);

	cells.push_back(p_cell);



	// Make cell population (2D)
	NodeBasedCellPopulation<2> cell_population(mesh, cells);
	cell_population.AddPopulationWriter<NodeLocationWriter>();
	//cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files

	// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
	double widthOfDomain = 10;
	// lengthOfDomain is a parameter
	ChastePoint<2> lower(-widthOfDomain*0.5, 0);
	ChastePoint<2> upper(widthOfDomain*0.5, lengthOfDomain);
	MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

	// Make PDE
	//UniformSourceEllipticPde<3> pde(0);
	MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (0.0));
	MAKE_PTR_ARGS(FunctionalBoundaryCondition<2>, p_functional_bc, (&bc_func));
	//MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
	bool is_neumann_bc = false; // Dirichlet BCs

	// Create a PDE modifier and set the name of the dependent variable in the PDE
	MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_functional_bc, is_neumann_bc, p_cuboid, 1.0));
	p_pde_modifier->SetDependentVariableName("csf1");
	p_pde_modifier->SetOutputGradient(true);
	p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

	/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
	 * (this time with dimension 2) and set the output directory, output multiple and end time. */
	OffLatticeSimulation<2> simulator(cell_population);
	simulator.AddSimulationModifier(p_pde_modifier);
	// Create and set an output directory that is different for each simulation
	std::stringstream output_directory;
	output_directory << "AddingMacrophages/ChemotaxisImplementationTesting/yLinearBCs/" << id_string;
	simulator.SetOutputDirectory(output_directory.str());
	simulator.SetSamplingTimestepMultiple(12);
	simulator.SetEndTime(120.0);


	/* Again we create a force law (this time with dimension 2), and pass it to the {{{OffLatticeSimulation}}}.*/
//	MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForceWithMacrophages<2>, p_force);
	// Anything within the set distance of a macrophage has half the spring force of a usual tumour-tumour spring.
//	p_force->SetHomotypicMacrophageSpringConstantMultiplier(1.0); // This is irrelvant in a test with only 1 macrophage
//	p_force->SetLabelledMacrophageSpringConstantMultiplier(labelledMacrophageSpringConstantMultiplier);
//	p_force->SetHomotypicLabelledSpringConstantMultiplier(homotypicLabelledSpringConstantMultiplier);
//	p_force->SetHeterotypicSpringConstantMultiplier(heterotypicSpringConstantMultiplier);
//
//	p_force->SetMeinekeSpringStiffness(meinekeSpringStiffness);
//	p_force->SetCutOffLength(1.5);

//	simulator.AddSimulationModifier(p_macProx);
//	simulator.AddForce(p_force);



	MAKE_PTR(DiffusionForce<2>, p_diffusion_force);
    /*
     * Compute the diffusion coefficient D as D = k*T/(6*pi*eta*r), where
     *
     * k = Boltzmann's constant,
     * T = absolute temperature,
     * eta = dynamic viscosity,
     * r = cell radius = 1
     */
	double D = 1.0; // Desired Diffusion Constant

	//p_diffusion_force->SetViscosity();
	p_diffusion_force->SetAbsoluteTemperature(temperatureKelvin);
	simulator.AddForce(p_diffusion_force);


	// Chemotaxis
	MAKE_PTR(ChemotacticForceCSF1<2>, p_chemotactic_force);
	p_chemotactic_force->SetChemotaxisSensitivity(chemotaxisSensitivity);
	simulator.AddForce(p_chemotactic_force);

	/* To run the simulation, we call {{{Solve()}}}. */

	try
	{
		simulator.Solve();

		OutputFileHandler results_handler(output_directory.str(), false);
		out_stream results_file = results_handler.OpenOutputFile("results.completiontimeandlocation");

		Node<2>* p_MacNode = cell_population.GetNode(nodeNum);
		const c_vector<double, 2>& r_MacNode_location = p_MacNode->rGetLocation();

		// Output summary statistics to results file
		(*results_file) << "ID" << ","
						<< "Simulation error" << ","
						<< "Macrophage x coordinate at end" << ","
						<< "Macrophage y coordinate at end" << ","
						<< "Simulation time at end" << std::endl;

		(*results_file) << id_string << ","
						<< "false" << ","
						<< boost::lexical_cast<std::string>(r_MacNode_location[0]) << ","
						<< boost::lexical_cast<std::string>(r_MacNode_location[1]) << ","
						<< boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTime());
		// Tidy up
		results_file->close();
	}
	catch(...)//std::string error)
	{
		//std::cout << error;
		OutputFileHandler results_handler(output_directory.str(), false);
		out_stream results_file = results_handler.OpenOutputFile("results.completiontimeandlocation");

		Node<2>* p_MacNode = cell_population.GetNode(nodeNum);
		const c_vector<double, 2>& r_MacNode_location = p_MacNode->rGetLocation();

		// Output summary statistics to results file
		(*results_file) << "ID" << ","
						<< "Simulation error" << ","
						<< "Macrophage x coordinate at end" << ","
						<< "Macrophage y coordinate at end" << ","
						<< "Simulation time at end" << std::endl;

		(*results_file) << id_string << ","
						<< "true" << ","
						<< boost::lexical_cast<std::string>(r_MacNode_location[0]) << ","
						<< boost::lexical_cast<std::string>(r_MacNode_location[1]) << ","
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
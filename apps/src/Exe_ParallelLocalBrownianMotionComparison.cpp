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
#include "ExecutableSupport.hpp"

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

#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

#include "LocalBrownianMotion.hpp"
#include "DiffusionForceChooseD.hpp"
#include "MacrophageProximityLabellingModifier.hpp"

#include "CuboidPeriodicBoundaryCondition.hpp"

#include "UniformSourceEllipticPde.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"

#include "ChemotacticForceCSF1.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



/*
 * Prototype functions
 */
void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(std::string id_string, double diffusionCoefficient, bool isLocalBM, double chemotaxisSensitivity, double radiusOfresolution, double iterationNumber);
void OutputToConsole(std::string id_string, std::string leading);
double bc_func(const ChastePoint<2>& p);


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
													("D", boost::program_options::value<double>()->default_value(1.0),"diffusionCoefficient: coefficient of diffusion")
													("Local", boost::program_options::value<bool>()->default_value(false),"isLocalBM: use local Brownian Motion or Global")
													("CS", boost::program_options::value<double>()->default_value(1.0),"chemotaxisSensitivity: Sensitivity of MP to chemotactic gradient")
													("RR", boost::program_options::value<double>()->default_value(5.0),"radiusOfResolution: Radius around each macrophage within which all cells are subject to Brownian motion")
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
	double diffusionCoefficient = variables_map["D"].as<double>();
	double isLocalBM = variables_map["Local"].as<bool>();
	double chemotaxisSensitivity = variables_map["CS"].as<double>();
	double radiusOfResolution = variables_map["RR"].as<double>();
	double iterationNumber = variables_map["IN"].as<double>();

	OutputToConsole(id_string, "Started");

	SetupSingletons();
	SetupAndRunSimulation(id_string,diffusionCoefficient,isLocalBM,chemotaxisSensitivity,radiusOfResolution,iterationNumber);
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
	return p[1]; // Return y value on boundary - yLinearBCs
}

void SetupAndRunSimulation(std::string id_string, double diffusionCoefficient, bool isLocalBM, double chemotaxisSensitivity, double radiusOfResolution, double iterationNumber)
{

	// Test macrophage on constant chemotactic gradient and log time to traverse.
	/* We vary parameters:
	 *
	 * diffusionCoefficient: D for Brownian motion
	 * isLocalBM: use local Brownian motion, or global?
	 * chemotaxisSensitivity: Sensitivity of MP to chemotactic gradient
	 * radiusOfResolution: radius around each macrophage within which all cells are subject to Brownian motion
	 * iterationNumber: e.g, simulation 4 of a given n for each parameter set (specified in python) - not used in simulation, just as an iterator
	 *
	 */
	// Generate Mesh:
	std::vector<Node<2>*> nodes;

	double lengthOfCuboid = 30;
	double widthOfCuboid = 30;

	unsigned nodeNum=0;
	for (double y=1; y<lengthOfCuboid; y++)
	{
		for (double x=1; x<widthOfCuboid; x++)
		{
			nodes.push_back(new Node<2>(nodeNum,  false,  x, y));
			nodeNum++;
		}
	}
	// Macrophage node
	//nodes.push_back(new Node<2>(nodeNum,  false,  widthOfCuboid*0.5, 0));
	// For central (non-gradient) test - add 0.1 to stop it being colocated with tumour cell above!
	nodes.push_back(new Node<2>(nodeNum,  false,  widthOfCuboid*0.5+0.1, lengthOfCuboid*0.5 + 0.1));

	NodesOnlyMesh<2> mesh;
	// Cut off length: 1.5 cell radii
	mesh.ConstructNodesWithoutMesh(nodes, 1.5);

	// Make cell pointers
	std::vector<CellPtr> cells;
	//MAKE_PTR(WildTypeCellMutationState, p_state);
	boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
	MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
	MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

	// Loop over nodes and create cells manually
	for (unsigned i=0; i<(nodes.size()-1); i++)
	{
		NoCellCycleModel* p_model = new NoCellCycleModel;
		p_model->SetDimension(2);

		CellPtr p_cell(new Cell(p_state, p_model));
		p_cell->SetCellProliferativeType(p_differentiated_type);

		cells.push_back(p_cell);
	}

	// Make Macrophage
	NoCellCycleModel* p_model = new NoCellCycleModel;
	p_model->SetDimension(2);

	CellPtr p_cell(new Cell(p_state, p_model));
	p_cell->SetCellProliferativeType(p_macrophage_type);

	cells.push_back(p_cell);

	// Make cell population (2D)
	NodeBasedCellPopulation<2> cell_population(mesh, cells);
	cell_population.AddPopulationWriter<NodeLocationWriter>();

	// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE (bigger than cells)
	ChastePoint<2> lower(0, -5);
	ChastePoint<2> upper(widthOfCuboid, lengthOfCuboid);
	MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

	// Make PDE
	MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (0.0));
	MAKE_PTR_ARGS(FunctionalBoundaryCondition<2>, p_functional_bc, (&bc_func));
	bool is_neumann_bc = false; // Dirichlet BCs

	// Create a PDE modifier and set the name of the dependent variable in the PDE
	MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_functional_bc, is_neumann_bc, p_cuboid, 1.0));
	p_pde_modifier->SetDependentVariableName("csf1");
	p_pde_modifier->SetOutputGradient(true);
	p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

	OffLatticeSimulation<2> simulator(cell_population);
	simulator.AddSimulationModifier(p_pde_modifier);

	// Create and set an output directory
	std::stringstream output_directory;
	output_directory << "CoarseBrownianMotion/MonolayerInfiltrationUnderCSFGradient_2D_6/" << id_string;
	simulator.SetOutputDirectory(output_directory.str());
	simulator.SetSamplingTimestepMultiple(60); // Was 120
	simulator.SetEndTime(50);

	// Add periodic boundary conditions, but away from top and bottom of cells
	MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<2>,p_boundary_condition,(&cell_population,lower,upper));
	simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);


	if(isLocalBM)
	{
		// Before adding force law, we need the helper modifier MacrophageProximityLabellingModifier to label cells close to macrophages
		MAKE_PTR(MacrophageProximityLabellingModifier<2>, p_macProx);

		MAKE_PTR(LocalBrownianMotion<2>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(diffusionCoefficient);
		// Specify radius within which cells are deemed close to a macrophage
		p_diffusion_force->SetRadiusOfResolution(radiusOfResolution);
		// Update labeller with correct radius
		p_macProx->SetMacrophageProximityLabellingRadius(radiusOfResolution);

		simulator.AddSimulationModifier(p_macProx);
		simulator.AddForce(p_diffusion_force);
	}
	if(not isLocalBM)
	{
		// Just use normal BM
		MAKE_PTR(DiffusionForceChooseD<2>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(diffusionCoefficient);

		simulator.AddForce(p_diffusion_force);
	}
	/* Create a force law, and pass it to the {{{OffLatticeSimulation}}}.*/
	MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
	p_force->SetMeinekeSpringStiffness(30.0);
	p_force->SetCutOffLength(1.5);
	simulator.AddForce(p_force);

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
						<< "Macrophage x coordinate at end" << ",";
		(*results_file) << "Macrophage y coordinate at end" << ",";

		(*results_file) << "Simulation time at end" << std::endl;

		(*results_file) << id_string << ","
						<< "false" << ","
						<< boost::lexical_cast<std::string>(r_MacNode_location[0]) << ",";
		(*results_file) << boost::lexical_cast<std::string>(r_MacNode_location[1]) << ",";

		(*results_file) << boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTime());
		// Tidy up
		results_file->close();
	}
	catch(Exception& e)//std::string error)
	{
		//throw(e);
		//std::cout << error;
		OutputFileHandler results_handler(output_directory.str(), false);
		out_stream results_file = results_handler.OpenOutputFile("results.completiontimeandlocation");

		Node<2>* p_MacNode = cell_population.GetNode(nodeNum);
		const c_vector<double, 2>& r_MacNode_location = p_MacNode->rGetLocation();

		// Output summary statistics to results file
		(*results_file) << "ID" << ","
						<< "Simulation error" << ","
						<< "Macrophage x coordinate at end" << ",";
		(*results_file) << "Macrophage y coordinate at end" << ",";

		(*results_file) << "Simulation time at end" << std::endl;

		(*results_file) << id_string << ","
						<< "true" << ","
						<< boost::lexical_cast<std::string>(r_MacNode_location[0]) << ",";
		(*results_file) << boost::lexical_cast<std::string>(r_MacNode_location[1]) << ",";

		(*results_file) << boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTime());
		// Tidy up
		results_file->close();
	}


	/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
	for (unsigned i=0; i<nodes.size(); i++)
	{
		delete nodes[i];
	}

}

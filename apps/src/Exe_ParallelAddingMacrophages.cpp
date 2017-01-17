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


/*
 * Prototype functions
 */
void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(std::string id_string, double csfProduction,double quiescentConcentration,double hypoxicConcentration, unsigned temperatureKelvin);
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
									("csfProd", boost::program_options::value<double>()->default_value(0.03),"csfProduction: The production rate of CSF by quiescent tumour cells")
									("QC", boost::program_options::value<double>()->default_value(0.7),"Quiescent Concentration: The threshold oxygen concentration for onset of quiescence")
									("HC", boost::program_options::value<double>()->default_value(0.5),"Hypoxic Concentration: The threshold oxygen concentration for onset of hypoxia")
									("TK", boost::program_options::value<unsigned>()->default_value(9),"Temperature in Kelvin: For the diffusion rate");



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
	double csfProduction = variables_map["csfProd"].as<double>();
	double quiescentConcentration = variables_map["QC"].as<double>();
	double hypoxicConcentration = variables_map["HC"].as<double>();
	double temperatureKelvin = variables_map["TK"].as<unsigned>();

	OutputToConsole(id_string, "Started");

	SetupSingletons();
	SetupAndRunSimulation(id_string,csfProduction,quiescentConcentration,hypoxicConcentration,temperatureKelvin);
	DestroySingletons();

	OutputToConsole(id_string, "Completed");
}

void SetupSingletons()
{
	// Set up what the test suite would do
	SimulationTime::Instance()->SetStartTime(0.0);

	// Reseed with 0 for same random numbers each time, or time(NULL) or simulation_id to change each realisation
	RandomNumberGenerator::Instance()->Reseed(0);
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

void SetupAndRunSimulation(std::string id_string, double csfProduction,double quiescentConcentration,double hypoxicConcentration, unsigned temperatureKelvin)
{
	// In previous tests, hypoxic cells continue their cell cycle until they reach the G1 phase. Here, we modify the cell cycle to cause them to freeze wherever they are in their cell cycle.
	/** The next line is needed because we cannot currently run node based simulations in parallel. */
	//EXIT_IF_PARALLEL;

	// Generate Mesh:
	// Make Vector
	std::vector<Node<3>*> nodes;
	// Add some nodes

	// Tumour Nodes
	unsigned nodeNum=0;
	double initialRadius = 8.0;
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
	unsigned firstMacrophageNode = nodeNum;


	// Macrophage Nodes - in a shell at fixed radius
	unsigned numMacrophages=100;
	double macrophageSphereRadius = 10.0;

	std::vector<double> xVec;
	std::vector<double> yVec;
	std::vector<double> zVec;
	unsigned i=0;
	double random_x;
	double random_y;
	double zCoord;
	while(i < numMacrophages)
	{
		// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
		random_x = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
		random_y = macrophageSphereRadius*(2*(RandomNumberGenerator::Instance()->ranf()-0.5));
		xVec.push_back(random_x);
		yVec.push_back(random_y);

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
		zVec.push_back(zCoord);
		i++;
	}

	for(unsigned i=0; i<xVec.size();i++)
	{
		nodes.push_back(new Node<3>(nodeNum,  false,  xVec[i], yVec[i], zVec[i]));
		nodeNum++;
	}

	NodesOnlyMesh<3> mesh;
	// Cut off length: 1.5 cell radii
	mesh.ConstructNodesWithoutMesh(nodes, 1.5);

	// Make cell pointers
	std::vector<CellPtr> cells;
	//MAKE_PTR(WildTypeCellMutationState, p_state);
	boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
	MAKE_PTR(StemCellProliferativeType, p_stem_type);
	MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);


	// Create tumour cells manually
	for (unsigned i=0; i<firstMacrophageNode; i++)
	{
		SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic* p_model = new SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic;
		p_model->SetDimension(3);
		//p_model->SetQuiescentConcentration(0.8); // Arbitrary
		//p_model->SetHypoxicConcentration(0.7); // Arbitrary

		CellPtr p_cell(new Cell(p_state, p_model));
		p_cell->SetCellProliferativeType(p_stem_type);
		p_cell->GetCellData()->SetItem("oxygen", 1.0);
		p_cell->GetCellData()->SetItem("csf1", 1.0);

		p_model->SetStemCellG1Duration(8.0);
		p_model->SetTransitCellG1Duration(8.0);

		p_model->SetStemCellG1Duration(4.0);
		p_model->SetHypoxicConcentration(hypoxicConcentration);
		p_model->SetQuiescentConcentration(quiescentConcentration);
		p_model->SetCriticalHypoxicDuration(8);//was 8

		double birth_time = - RandomNumberGenerator::Instance()->ranf() *
				(  p_model->GetStemCellG1Duration()
						+ p_model->GetSG2MDuration() );
		p_cell->SetBirthTime(birth_time);
		cells.push_back(p_cell);
	}
	// Now create macrophages
	for (unsigned i=firstMacrophageNode; i<nodes.size(); i++)
	{
		// Make Macrophage
		NoCellCycleModel* p_model2 = new NoCellCycleModel;
		p_model2->SetDimension(3);

		CellPtr p_cell(new Cell(p_state, p_model2));
		p_cell->SetCellProliferativeType(p_macrophage_type);
		p_cell->GetCellData()->SetItem("csf1", 1.0);

		cells.push_back(p_cell);
	}

	// Make cell population (3D)
	NodeBasedCellPopulation<3> cell_population(mesh, cells);
	cell_population.AddPopulationWriter<VoronoiDataWriter>(); // Write VTU files
	// Write summary statistic files
	cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
	cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
	cell_population.AddCellWriter<CellMutationStatesWriter>();
	cell_population.AddCellWriter<CellProliferativeTypesWriter>();

	// Make PDE (Oxygen)
	MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.03));//was -0.03
	MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
	bool is_neumann_bc = false; // Dirichlet BCs

	// Make PDE (CSF1)
	MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pdeCSF, (cell_population, csfProduction));
	//MAKE_PTR_ARGS(AveragedSourceParabolicPde<3>, p_pdeCSF, (cell_population, 1.0, 1.0, 0.03));//du/dt coefficient, diffusion coefficient, source coefficient
	MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bcCSF, (0.0));
	bool is_neumann_bcCSF = true; // Neumann BCs


	// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
	double cubeDomainDistanceToBoundary = 20;
	c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
	ChastePoint<3> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
	ChastePoint<3> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
	MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

	// Create a PDE modifier and set the name of the dependent variable in the PDE
	MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
	//MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<3>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 2.0));
	p_pde_modifier->SetDependentVariableName("oxygen");
	p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

	MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifierCSF, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 2.0));
	p_pde_modifierCSF->SetDependentVariableName("csf1");
	p_pde_modifierCSF->SetOutputGradient(true);
	p_pde_modifierCSF->SetBcsOnBoxBoundary(true); //was false




	/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
	 * (this time with dimension 3) and set the output directory, output multiple and end time. */
	OffLatticeSimulation<3> simulator(cell_population);
	simulator.AddSimulationModifier(p_pde_modifier);
	simulator.AddSimulationModifier(p_pde_modifierCSF);


	// Create and set an output directory that is different for each simulation
	std::stringstream output_directory;
	output_directory << "AddingMacrophages/InfiltrationIntoSpheroidParameterSweep/" << id_string;
	simulator.SetOutputDirectory(output_directory.str());

	simulator.SetSamplingTimestepMultiple(12);//12);
	simulator.SetEndTime(30.0);

	/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
	MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
	p_force->SetMeinekeSpringStiffness(30.0);
	p_force->SetCutOffLength(1.5);
	simulator.AddForce(p_force);

	MAKE_PTR(DiffusionForce<3>, p_diffusion_force);
	//p_diffusion_force->SetViscosity();
	p_diffusion_force->SetAbsoluteTemperature(temperatureKelvin);
	simulator.AddForce(p_diffusion_force);

	// Add Chemotaxis
	MAKE_PTR(ChemotacticForceCSF1<3>, p_chemotactic_force);
	simulator.AddForce(p_chemotactic_force);

	//// Now put in a random killer
	//double probOfDeathPerHour;
	//probOfDeathPerHour = 0.01;
	//MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer, (&cell_population, probOfDeathPerHour));
	//simulator.AddCellKiller(p_cell_killer);

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

#ifndef TESTDORIE1982REPRODUCERESULTS_HPP_
#define TESTDORIE1982REPRODUCERESULTS_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

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
//#include "Debug.hpp"

#include "AveragedSourceEllipticPde.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
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

#include "CuboidPeriodicBoundaryCondition.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



class TestDorie1982ReproduceResults : public AbstractCellBasedTestSuite
{
public:

	void TestSpheroidGrowthCreateBenchmarkTumour() throw(Exception)
    															{
		/* This is a simulation designed to provide groundwork to reproducing Dorie et al.
		 * "Migration and Internalization of Cells and Polystyrene Microspheres in Tumour Cell Spheroids"
		 * Experimental Cell Research 141 (1982) 201-209
		 *
		 * In this case, inert polystyrene beads infiltrate a tumour spheroid where the only forces at work are proliferation and brownian motion
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<3>*> nodes;
		// Add some nodes

		double cubeDomainDistanceToBoundary = 25;

		// Tumour Nodes
		unsigned nodeNum=0;
		double initialRadius = 2;//13.5;
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

		NodesOnlyMesh<3> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

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

			p_model->SetHypoxicConcentration(0.15);
			p_model->SetQuiescentConcentration(0.35);
			p_model->SetCriticalHypoxicDuration(8);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(  p_model->GetStemCellG1Duration()
							+ p_model->GetSG2MDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}


		// Make cell population (3D)
		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
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
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("oxygen");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		simulator.SetOutputDirectory("FullSimulations/Dorie1982/SpheroidGrowth/IronManTest");
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(1.0);

		// Add periodic boundary conditions for cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<3>,p_periodic_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_periodic_boundary_condition);
		//TODO Archiving Problem

		// Add Brownian motion for all cells
		MAKE_PTR(DiffusionForceChooseD<3>, p_diffusion_force);
//		MAKE_PTR(DiffusionForce<3>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.05);
		simulator.AddForce(p_diffusion_force);
		//TODO Archiving Problem in "choose D"

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

		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
    															}

	void TestSpheroidGrowthLoadBenchmarkTumour() throw(Exception)
	    															{
		/* This is a simulation designed to provide groundwork to reproducing Dorie et al.
		 * "Migration and Internalization of Cells and Polystyrene Microspheres in Tumour Cell Spheroids"
		 * Experimental Cell Research 141 (1982) 201-209
		 *
		 * In this case, inert polystyrene beads infiltrate a tumour spheroid where the only forces at work are proliferation and brownian motion
		 */

		OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("FullSimulations/Dorie1982/SpheroidGrowth/IronManTest", 1.0);
		p_simulator->SetEndTime(2.0);

		// We remake the cuboid for the PDE FE Mesh
		double cubeDomainDistanceToBoundary = 25;
		ChastePoint<3> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary);
		ChastePoint<3> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

		// Now loop over simulation modifiers
        for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
             iter != p_simulator->GetSimulationModifiers()->end();
             ++iter)
        {

        	// Attempt to convert modifier to EllipticBoxDomainPdeModifier. If modifier is of a different type, this should return NULL
        	boost::shared_ptr<EllipticBoxDomainPdeModifier<3> > p_newModifier = boost::dynamic_pointer_cast<EllipticBoxDomainPdeModifier<3> >(*iter);

        	// If not NULL, we have found the correct modifier
        	if(p_newModifier)
        	{
        		// Now regenerate the FeMesh.
        		p_newModifier->GenerateFeMesh(p_cuboid,1.0);
        	}
        }

		p_simulator->Solve();
		CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);
		delete p_simulator;



	    															}

};


#endif /* TESTDORIE1982REPRODUCERESULTS_HPP_ */

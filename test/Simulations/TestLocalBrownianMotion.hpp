#ifndef TESTLOCALBROWNIANMOTION_HPP_
#define TESTLOCALBROWNIANMOTION_HPP_


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

/**
 * Functional Boundary Condition, setting BC = bcConc if z = distance
 */
double bc_func(const ChastePoint<2>& p)
{
	return p[1];
}

class TestLocalBrownianMotion : public AbstractCellBasedTestSuite
{
public:

	void TestLocalBrownianMotionInCuboid2D() throw(Exception)
	{
		EXIT_IF_PARALLEL;

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
		nodes.push_back(new Node<2>(nodeNum,  false,  widthOfCuboid*0.5, 0));

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
		ChastePoint<2> lower(0, 0);
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
		output_directory << "CoarseBrownianMotion/Cuboid2D_CSFGradient/1";
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(100);

		// Add periodic boundary conditions, but away from top and bottom of cells
		MAKE_PTR_ARGS(CuboidPeriodicBoundaryCondition<2>,p_boundary_condition,(&cell_population,lower,upper));
		simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

		// Before adding force law, we need the helper modifier MacrophageProximityLabellingModifier to label cells close to macrophages
		MAKE_PTR(MacrophageProximityLabellingModifier<2>, p_macProx);

		MAKE_PTR(LocalBrownianMotion<2>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.1);
		// Specify radius within which cells are deemed close to a macrophage
		double radiusOfResolution = 10.0;
		p_diffusion_force->SetRadiusOfResolution(radiusOfResolution);
		// Update labeller with correct radius
		p_macProx->SetMacrophageProximityLabellingRadius(radiusOfResolution);

		simulator.AddSimulationModifier(p_macProx);
		simulator.AddForce(p_diffusion_force);

		/* Create a force law, and pass it to the {{{OffLatticeSimulation}}}.*/
		MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
		p_force->SetMeinekeSpringStiffness(30.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Chemotaxis
		MAKE_PTR(ChemotacticForceCSF1<2>, p_chemotactic_force);
		p_chemotactic_force->SetChemotaxisSensitivity(2.0);
		simulator.AddForce(p_chemotactic_force);

		simulator.Solve();


	}

};


#endif /* TESTLOCALBROWNIANMOTION_HPP_ */

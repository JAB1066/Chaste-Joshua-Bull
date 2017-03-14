#ifndef TESTMACROPHAGECHEMOTAXIS_HPP_
#define TESTMACROPHAGECHEMOTAXIS_HPP_


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
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines a node-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"
// Write VTU files
#include "VoronoiDataWriter.hpp"
#include "RandomCellKiller.hpp"

// For test 3
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

#include "NodeLocationWriter.hpp"
#include "ChemotacticForceCSF1.hpp"
#include "DiffusionForceChooseD.hpp"

/**
 * Functional Boundary Condition, setting BC = bcConc if z = distance
 */
double bc_func(const ChastePoint<3>& p)
{
	double bcConc = 1.0;
	double distance = 10;
	double value = 0.0;

	if(p[2] == distance)
	{
		value = bcConc;
	}

	return value;
}

double bc_func2(const ChastePoint<2>& p)
{
	double bcConc = 1.0;
	double distance = 15;
	double value = 0.0;

	if(p[1] == distance)
	{
		value = bcConc;
	}

	return value;
}

class TestMacrophageChemotaxis : public AbstractCellBasedTestSuite
{
public:

	void TestMacrophageFollowingTumourCell() throw(Exception)
													{
		// This test simulates a moving tumour cell generating a CSF-1 gradient, which attracts a macrophage.
		EXIT_IF_PARALLEL;

		// Make Vector
		std::vector<Node<2>*> nodes;
		// Add macrophage node and TC node

		nodes.push_back(new Node<2>(0,  false,  10, 10));
		nodes.push_back(new Node<2>(1,  false,  -10, -10));

		NodesOnlyMesh<2> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;

		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);

		//boost::shared_ptr<AbstractCellProperty> p_prolifType = CellPropertyRegistry::Instance()->Get<MacrophageCellProliferativeType>();

		// Create cells manually
		// Macrophage:
		NoCellCycleModel* p_model = new NoCellCycleModel;
		p_model->SetDimension(2);
		CellPtr p_cell(new Cell(p_state, p_model));
		p_cell->SetCellProliferativeType(p_macrophage_type);
		p_cell->GetCellData()->SetItem("csf1", 1.0);
		cells.push_back(p_cell);

		// Tumour cell
		p_model->SetDimension(2);
		p_cell->SetCellProliferativeType(p_stem_type);
		p_cell->GetCellData()->SetItem("csf1", 1.0);
		cells.push_back(p_cell);

		// Make cell population (2D)
		NodeBasedCellPopulation<2> cell_population(mesh, cells);
		cell_population.AddPopulationWriter<NodeLocationWriter>(); // Write VTU files
		// Write summary statistic files
		cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
		cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		double squareDomainDistanceToBoundary = 75;
		//c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
		c_vector<double,2> centroid;
		centroid(0)=0;
		centroid(1)=0;
		ChastePoint<2> lower(centroid(0)-squareDomainDistanceToBoundary, centroid(1)-squareDomainDistanceToBoundary);
		ChastePoint<2> upper(centroid(0)+squareDomainDistanceToBoundary, centroid(1)+squareDomainDistanceToBoundary);
		MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

		// Make PDE (CSF1)
		//MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pdeCSF, (cell_population, 30.0, 1.0));
		MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pdeCSF, (cell_population, 1.0, 1.0, 5.0));//du/dt coefficient, diffusion coefficient, source coefficient
		MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bcCSF, (0.0));
		bool is_neumann_bcCSF = false; // Dirichlet BCs


		// Create a PDE modifier and set the name of the dependent variable in the PDE
		//MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 2.0));
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pdeCSF, p_bcCSF, is_neumann_bcCSF, p_cuboid, 1.0));
		p_pde_modifier->SetDependentVariableName("csf1");
		p_pde_modifier->SetOutputGradient(true);
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<2> simulator(cell_population);
		simulator.AddSimulationModifier(p_pde_modifier);
		//simulator.AddForce(p_diffusion_force);
		simulator.SetOutputDirectory("AddingMacrophages/MacrophageFollowingTC");
		simulator.SetSamplingTimestepMultiple(12);//12);
		simulator.SetEndTime(2000.0);


		// Make Chemotactic Force for CSF1
		MAKE_PTR(ChemotacticForceCSF1<2>, p_chemotactic_force);
		p_chemotactic_force->SetChemotaxisSensitivity(1.0);
		simulator.AddForce(p_chemotactic_force);

		// Make diffusive force
		MAKE_PTR(DiffusionForceChooseD<2>, p_diffusion_force);
		p_diffusion_force->SetDiffusionScalingConstant(0.05);
		simulator.AddForce(p_diffusion_force);


		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
													}

};


#endif /* TESTMACROPHAGECHEMOTAXIS_HPP_ */

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "HeterotypicBoundaryLengthWriter.hpp"

#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
//#include "RandomMotionForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
//#include "RandomMotionForce.hpp"

#include "OnLatticeSimulation.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"

#include "CaBasedCellPopulation.hpp"
#include "ShovingCaBasedDivisionRule.hpp"
#include "RandomCaSwitchingUpdateRule.hpp"
#include "DifferentialAdhesionCaSwitchingUpdateRule.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"

#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

static const double M_TIME_TO_STEADY_STATE = 10; //10
static const double M_TIME_FOR_SIMULATION = 100; //100
static const double M_NUM_CELLS_ACROSS = 20; //20 // this ^2 cells
static const double M_CELL_FLUCTUATION = 1.0;

class TestCellSortingLiteratePaper : public AbstractCellBasedWithTimingsTestSuite
{
private:
	void RandomlyLabelCells(std::list<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
	    {
	        for (std::list<CellPtr>::iterator cell_iter = rCells.begin();
	             cell_iter != rCells.end();
	             ++cell_iter)
	        {
	            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
	            {
	               (*cell_iter)->AddCellProperty(pLabel);
	            }
	        }
	    }

	public:

	void TestCaBasedMonolayerCellSorting() throw (Exception)
	    {
	        // Create a simple 2D PottsMesh
	        unsigned domain_wide = 2*M_NUM_CELLS_ACROSS;

	        PottsMeshGenerator<2> generator(domain_wide, 0, 0, domain_wide, 0, 0);
	        PottsMesh<2>* p_mesh = generator.GetMesh();

	        p_mesh->Translate(-(double)domain_wide*0.5 + 0.5,-(double)domain_wide*0.5 + 0.5);

	        // Specify where cells lie
	        std::vector<unsigned> location_indices;
	        for (unsigned i=0; i<M_NUM_CELLS_ACROSS; i++)
	        {
	          for (unsigned j=0; j<M_NUM_CELLS_ACROSS; j++)
	          {
	              unsigned offset = (domain_wide+1) * (domain_wide-M_NUM_CELLS_ACROSS)/2;
	              location_indices.push_back(offset + j + i * domain_wide );
	          }
	        }
	        std::vector<CellPtr> cells;
	        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
	        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
	        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

	        // Create cell population
	        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

	        // Set population to output all data to results files
	        cell_population.AddCellWriter<CellIdWriter>();
	        cell_population.AddCellWriter<CellMutationStatesWriter>();
	        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
	        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

	        OnLatticeSimulation<2> simulator(cell_population);
	        simulator.SetOutputDirectory("CellSorting/Ca");
	        simulator.SetDt(0.01);
	        simulator.SetSamplingTimestepMultiple(100);
	        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

	        // Add Division Rule
	        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule(new ShovingCaBasedDivisionRule<2>());
	        cell_population.SetCaBasedDivisionRule(p_division_rule);

	        // Add switching Update Rule
	        MAKE_PTR(DifferentialAdhesionCaSwitchingUpdateRule<2u>, p_switching_update_rule);
	        p_switching_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.1);
	        p_switching_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.2);
	        p_switching_update_rule->SetCellCellAdhesionEnergyParameter(0.1);
	        p_switching_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.2);
	        p_switching_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.4);
	        p_switching_update_rule->SetTemperature(0.1);
	        simulator.AddUpdateRule(p_switching_update_rule);

	        simulator.Solve();

	        // Now label some cells
	        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
	        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

	        // modify parameters
	        p_switching_update_rule->SetTemperature(0.1*M_CELL_FLUCTUATION);

	        // Run simulation
	        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
	        simulator.Solve();

	        // Check that the same number of cells
	        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

	        // Test no births or deaths
	        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
	        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
	    }
}





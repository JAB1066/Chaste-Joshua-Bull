#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <iomanip>
#include <boost/foreach.hpp>
#include "OffLatticeSimulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
//#include "AveragedSourcePde.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "ConstBoundaryCondition.hpp"
//#include "CellBasedPdeHandler.hpp"
#include "ChastePoint.hpp"
#include "SmartPointers.hpp"
#include "ApoptoticCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestSpheroidExperimentsLiteratePaper : public AbstractCellBasedTestSuite
{
private:

    void setUp()
    {
        AbstractCellBasedTestSuite::setUp();
        CellBasedEventHandler::Reset();
    }
    void tearDown()
    {
        AbstractCellBasedTestSuite::tearDown();

        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();
    }

public:

    void TestMeshBasedSpheroidWithPde() throw(Exception)
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        //nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        // We use the stem cell G1 duration, so make these 'stem' cells
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            StochasticOxygenBasedCellCycleModel* p_model = new StochasticOxygenBasedCellCycleModel();
            p_model->SetDimension(3);
            p_model->SetStemCellG1Duration(4.0);
            p_model->SetHypoxicConcentration(0.1);
            p_model->SetQuiescentConcentration(0.3);
            p_model->SetCriticalHypoxicDuration(8);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (  p_model->GetStemCellG1Duration()
                                 + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(birth_time);
            p_cell->SetCellProliferativeType(p_stem_type);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);
        cell_population.SetWriteVtkAsPoints(true);

        cell_population.SetDataOnAllCells("oxygen", 1.0);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetEndTime(1000); // hours

        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetOutputDirectory("Plos2013_MeshBasedSpheroidWithPde");


        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<3>, p_pde, (cell_population, -1.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));
        bool is_neumann_bc = false;

		MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc));
		p_pde_modifier->SetDependentVariableName("oxygen");
		simulator.AddSimulationModifier(p_pde_modifier);

        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetMeinekeSpringStiffness(30.0); // default is 15.0;
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(ApoptoticCellKiller<3>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);
    }

    void TestLongerMeshBasedSpheroidWithPde() throw(Exception)
    {
        FileFinder test_data_directory("Plos2013_MeshBasedSpheroidWithPde/archive",
                                       RelativeTo::ChasteTestOutput);

        OutputFileHandler archive_handler("Plos2013_LongerMeshBasedSpheroidWithPde/archive");

        // Following is done in two lines to avoid a bug in Intel compiler v12.0!
        std::vector<FileFinder> temp_files = test_data_directory.FindMatches("*");
        BOOST_FOREACH(FileFinder temp_file, temp_files)
        {
            archive_handler.CopyFileTo(temp_file);
        }

        OffLatticeSimulation<3>* p_simulator
            = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("Plos2013_LongerMeshBasedSpheroidWithPde", 1000);

        p_simulator->SetEndTime(1500);
        p_simulator->SetOutputDirectory("Plos2013_LongerMeshBasedSpheroidWithPde");

        p_simulator->Solve();

        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);

    }
};

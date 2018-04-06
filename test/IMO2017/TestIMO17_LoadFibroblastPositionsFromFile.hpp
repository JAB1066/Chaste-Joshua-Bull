#ifndef TEST_IMO17_FIBROBLASTTUMOURPATTERNS_HPP_
#define TEST_IMO17_FIBROBLASTTUMOURPATTERNS_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
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

//#include "SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
//#include "ChemotacticForceCSF1.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

#include "AveragedSourceParabolicPde.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ApoptoticCellKiller.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "CuboidPeriodicBoundaryCondition.hpp"
#include "VolumeTrackingModifier.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "RandomCellKiller.hpp"
#include "AbstractCellProperty.hpp"

#include "MacrophageProximityLabellingModifier.hpp"
#include "AddProliferationInhibitingTreatmentModifier.hpp"
#include "UniformContactInhibitionCellCycleModelSlowByFibroblasts.hpp"
#include "ResistantIMOMutationState.hpp"
#include "FibroblastIMOMutationState.hpp"
//#include "ResistantMutationStateTrackingModifier.hpp"
#include "RandomCellKillerIgnoringFibroblastsIMO.hpp"
#include "KeepMacrophagesStationaryLinearSpringForce.hpp"
#include "RepulsionForceKeepMacrophagesStationary.hpp"

#include "ColumnDataReader.hpp"



#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



class TestIMO17_LoadFibroblastPositionsFromFile : public AbstractCellBasedTestSuite
{
public:

    void TestAddingTreatmentRandomFibroblasts() throw(Exception)
	{
		/* We throw down fibroblasts randomly across the domain, and seed
		 * it with a tumour cell population. Fibroblasts produce growth factor GF, which
		 * promotes the growth of tumour cells.
		 *
		 * After time T, we simulate treatment: a "resistant" TC population is unaffected.
		 * The "non-resistant" cell population only proliferates in the presence of GF.
		 *
		 */

		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		// Parameters
		unsigned numberOfFibroblasts = 0;//150;
		unsigned numberOfSensitiveTumourCells = 250;
		unsigned numberOfResistantTumourCells = 25;
		double squareDomainWidth = 30;
		const int DIM = 2;
		double simulationDuration = 120;
		double timeToApplyTreatment = 10;
		double distanceFromFibroblastToFeelBenefit = 3.0;

		/** Treatment Parameters */
		// Pre-Treatment
		double MinCellCycleDuration_Sensitive = 0.1;
		double MaxCellCycleDuration_Sensitive = 0.2;

        double MinCellCycleDuration_Resistant = 0.2;
		double MaxCellCycleDuration_Resistant = 0.3;
		// Post-treatment
		double MinCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MaxCellCycleDuration_SensitivePostTreatment = DBL_MAX;
		double MinCellCycleDurationByFibroblast_SensitivePostTreatment = 0.1;
		double MaxCellCycleDurationByFibroblast_SensitivePostTreatment = 0.2;

		std::string outputDirectory = "IMO2017/SetDistributions/NoFibroblasts/";

		/*
        ifstream inFile;
        inFile.open("/mi/share/scratch/bull/ChasteStuff/chaste_test_output/IMO2017/INPUT_FibroblastDistributions/NoFibroblasts.csv");
        if (!inFile) {
            std::cout << "Unable to open file";
            exit(1); // terminate with error
        }

        std::string value;
        std::vector<std::string> record;
        while (inFile)
        {
            //std::cout << x << std::endl;
            std::getline(inFile, value,'\n');
            //std::cout << value << std::endl;

            std::istringstream ss( value );

            while (ss)
            {
                std::string s;
                if (!getline( ss, s, ',' )) break;
                record.push_back( s );
            }
            //sum = sum + x;
        }
        //std::cout << "Length of record is " << record.size() << std::endl;
        inFile.close();
        std::vector<double> xs;
        std::vector<double> ys;
        int index = 0;
        for(unsigned i = 0; i < numberOfFibroblasts*2; i=i+2)
        {
            //std::cout << "i = " << i << ", i/2 = " << index << std::endl;
            //std::cout << "Size of xs is " << xs.size() << std::endl;
            xs.push_back(stod(record[index]));
            index++;
            ys.push_back(stod(record[index]));
            index++;
        }
        */
		//double squareDomainWidthHalved = squareDomainWidth*0.5;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		// We randomly distribute nodes across the domain
		unsigned nodeNum=0;
		double x;
		double y;
		//double z;

		while(nodeNum < (numberOfSensitiveTumourCells + numberOfResistantTumourCells))
		{
			// Generate random x, y coordinates between -macrophageSphereRadius and + macrophageSphereRadius
			x = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			y = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf());//-squareDomainWidthHalved);
			//z = squareDomainWidth*(RandomNumberGenerator::Instance()->ranf()-squareDomainWidthHalved);

			nodes.push_back(new Node<DIM>(nodeNum,  false,  x, y));//, z));
			nodeNum++;
		}
		/*
		for(unsigned i = 0; i < numberOfFibroblasts; i++)
        {
            double x = xs[i];
            double y = ys[i];
            double scaledX = x*squareDomainWidth*1.0;
            double scaledY = y*squareDomainWidth*1.0;

            assert(scaledX < squareDomainWidth);
            assert(scaledY < squareDomainWidth);
            assert(scaledX > 0);
            assert(scaledY > 0);

            std::cout<<i<<" - x is "<<x<<", scaledX is "<<scaledX<<std::endl;
            std::cout<<i<<" - y is "<<y<<", scaledY is "<<scaledY<<std::endl;

            nodes.push_back(new Node<DIM>(nodeNum,  false,  scaledX, scaledY));//, z));
            nodeNum++;
        }
        double maxEle = *max_element(xs.begin(), xs.end());
        std::cout<<"Max value for xs : "<<maxEle<<std::endl;
        maxEle = *max_element(ys.begin(), ys.end());
        std::cout<<"Max value for ys : "<<maxEle<<std::endl;
*/
		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_wild_state);
		MAKE_PTR(ResistantIMOMutationState, p_resistant_state);
		MAKE_PTR(FibroblastIMOMutationState, p_fibroblast_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(MacrophageCellProliferativeType, p_fibroblast_type);

		// Create tumour cells manually
		for (unsigned i=0; i<numberOfSensitiveTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Sensitive);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Sensitive);

			CellPtr p_cell(new Cell(p_wild_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					( p_model->GetMinCellCycleDuration() );
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}
        for (unsigned i=0; i<numberOfResistantTumourCells; i++)
		{
			UniformContactInhibitionCellCycleModelSlowByFibroblasts* p_model = new UniformContactInhibitionCellCycleModelSlowByFibroblasts;
			p_model->SetDimension(DIM);
			p_model->SetQuiescentVolumeFraction(0.6);
			p_model->SetEquilibriumVolume(1.0);

			p_model->SetMinCellCycleDuration(MinCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDuration(MaxCellCycleDuration_Resistant);
			p_model->SetMinCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);
			p_model->SetMaxCellCycleDurationByFibroblast(MaxCellCycleDuration_Resistant);

			CellPtr p_cell(new Cell(p_resistant_state, p_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_model->GetMinCellCycleDuration();
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}

		for (unsigned i = 0; i<numberOfFibroblasts; i++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;

            CellPtr p_cell(new Cell(p_fibroblast_state, p_model));
            p_cell->SetCellProliferativeType(p_fibroblast_type);

            cells.push_back(p_cell);
        }


		// Make cell population (3D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.AddCellWriter<CellMutationStatesWriter>();
		//cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		ChastePoint<DIM> lower(0,0);
		ChastePoint<DIM> upper(squareDomainWidth,squareDomainWidth);

		//ChastePoint<3> lower(-squareDomainWidthHalved, -squareDomainWidthHalved);//, -squareDomainWidthHalved-0.1);
		//ChastePoint<3> upper(squareDomainWidthHalved, squareDomainWidthHalved);//, squareDomainWidthHalved+0.1);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));

		/*
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifierVariable<3>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc, p_cuboid, 2.0));
		p_pde_modifier->SetDependentVariableName("HGF");
		p_pde_modifier->SetBcsOnBoxBoundary(true); //was false
        */

		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
		 * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		//simulator.AddSimulationModifier(p_pde_modifier);

		simulator.SetOutputDirectory(outputDirectory);
		simulator.SetSamplingTimestepMultiple(24); // 12 timesteps = 6 minutes = 0.1 hours
		simulator.SetEndTime(simulationDuration);



		// Add boundary conditions
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = unit_vector<double>(2,1);
        normal(0) = -1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc1, (&cell_population, point, normal));
        p_bc1->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        normal(0) = 0.0; normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc2, (&cell_population, point, normal));
        p_bc2->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = squareDomainWidth; point(1) = squareDomainWidth;
        normal(0) = 0.0; normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc3, (&cell_population, point, normal));
        p_bc3->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        normal(0) = 1.0; normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<DIM>, p_bc4, (&cell_population, point, normal));
        p_bc4->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(AddProliferationInhibitingTreatmentModifier<DIM>, p_treatment_modifier);
        p_treatment_modifier->SetTreatmentStartTime(timeToApplyTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MinCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_PostTreatment(MaxCellCycleDuration_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MinCellCycleDurationByFibroblast_SensitivePostTreatment);
        p_treatment_modifier->SetMinCellCycleDuration_ByFibroblastPostTreatment(MaxCellCycleDurationByFibroblast_SensitivePostTreatment);
        simulator.AddSimulationModifier(p_treatment_modifier);

        // Add a cell killer
        double probabiltyOfDeathInOneHour = 0.05;
        MAKE_PTR_ARGS(RandomCellKillerIgnoringFibroblastsIMO<DIM>, p_random_killer, (&cell_population, probabiltyOfDeathInOneHour));
        simulator.AddCellKiller(p_random_killer);


		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		/* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
		//MAKE_PTR(GeneralisedLinearSpringForce<DIM>, p_force);
		//p_force->SetMeinekeSpringStiffness(15.0);
		//p_force->SetCutOffLength(1.5);
		//simulator.AddForce(p_force);

		//OR


        MAKE_PTR(RepulsionForceKeepMacrophagesStationary<DIM>, p_force);
        simulator.AddForce(p_force);

        /*
        MAKE_PTR(KeepMacrophagesStationaryLinearSpringForce<DIM>, p_force);
        p_force->SetMeinekeSpringStiffness(15.0);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);
        */




		// Label cells which are in close proximity to fibroblasts
		MAKE_PTR(MacrophageProximityLabellingModifier<DIM>, p_proximity_labeller);
		p_proximity_labeller->SetMacrophageProximityLabellingRadius(distanceFromFibroblastToFeelBenefit);
		simulator.AddSimulationModifier(p_proximity_labeller);

        /*** To run the simulation, we call {{{Solve()}}}. ***/
		simulator.Solve();

		//CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);


		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void dontTestReadFibroblastsFromData() throw(Exception)
	{
	    //int sum = 0;
        ifstream inFile;

        std::cout << "BKJDBAGKJBSG" << std::endl;
        inFile.open("/mi/share/scratch/bull/ChasteStuff/chaste_test_output/IMO2017/INPUT_FibroblastDistributions/CirclesFibs.csv");
        if (!inFile) {
            std::cout << "Unable to open file";
            exit(1); // terminate with error
        }

        std::string value;
        while (inFile)
        {
            //std::cout << x << std::endl;
            std::getline(inFile, value,'\n');
            std::cout << value << std::endl;

            std::istringstream ss( value );
            std::vector <std::string> record;

            while (ss)
            {
                std::string s;
                if (!getline( ss, s, ',' )) break;
                record.push_back( s );
            }
            //sum = sum + x;
        }

        inFile.close();
        //cout << "Sum = " << sum << endl;






        /*
		ifstream input_stream("IMO2017/INPUT_FibroblastDistributions/GridFibs.info");
        std::string value;
        std::cout << "BLAHAOSDGOSDNGOSDNGLSDKGNLAKN" << std::endl;
        bool x = true;
        while ((!input_stream.eof()) && x)
        {
            std::cout << "dsgsdngsdjgnsbmclblxn" << std::endl;
            std::getline(input_stream, value);
            std::cout << value;
            x = false;
        }
        */



	}
};



#endif /* TEST_IMO17_FIBROBLASTTUMOURPATTERNS_HPP_ */

/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "OutputOnlyMacrophageSummaryStatisticsModifier.hpp"
#include "OutputFileHandler.hpp"
#include "Debug.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "AbstractCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellData.hpp"

template<unsigned DIM>
OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::OutputOnlyMacrophageSummaryStatisticsModifier()
: AbstractCellBasedSimulationModifier<DIM>(),
  mOutputDirectory(""),
  mOutputFrequencyInHours(0.1),
  mQuiescenceLevel(0.0),
  mHypoxiaLevel(0.0),
  mTimestepsElapsed(0)
  {
  }

template<unsigned DIM>
OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::~OutputOnlyMacrophageSummaryStatisticsModifier()
{
}

template<unsigned DIM>
std::string OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::GetOutputDirectory()
{
	return mOutputDirectory;
}

template<unsigned DIM>
void OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
}

template<unsigned DIM>
double OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::GetOutputFrequencyInHours()
{
	return mOutputFrequencyInHours;
}

template<unsigned DIM>
void OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::SetOutputFrequencyInHours(double outputFrequencyInHours)
{
	mOutputFrequencyInHours = outputFrequencyInHours;
}

template<unsigned DIM>
double OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::GetQuiescenceLevel()
{
	return mQuiescenceLevel;
}

template<unsigned DIM>
void OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::SetQuiescenceLevel(double quiescenceLevel)
{
	mQuiescenceLevel = quiescenceLevel;
}

template<unsigned DIM>
double OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::GetHypoxiaLevel()
{
	return mHypoxiaLevel;
}

template<unsigned DIM>
void OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::SetHypoxiaLevel(double hypoxiaLevel)
{
	mHypoxiaLevel = hypoxiaLevel;
}

template<unsigned DIM>
void OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::WriteData(AbstractCellPopulation<DIM,DIM>& rCellPopulation, out_stream pOutStreamMacrophageNodes, out_stream pOutStreamCompartmentBoundaries, out_stream pOutStreamMacrophageOxygenConc, out_stream pOutStreamTumourCellCount)
{

	*pOutStreamMacrophageNodes << SimulationTime::Instance()->GetTime() << "\t";
	*pOutStreamCompartmentBoundaries << SimulationTime::Instance()->GetTime() << "\t";
	*pOutStreamMacrophageOxygenConc << SimulationTime::Instance()->GetTime() << "\t";
	*pOutStreamTumourCellCount << SimulationTime::Instance()->GetTime() << "\t";

	const c_vector<double,DIM>& centroid = rCellPopulation.GetCentroidOfCellPopulation();
	double HypoxiaBoundary = 0;
	double QuiescenceBoundary = 0;
	double TumourBoundary = 0;
	double InnerBoundary = 0;
	double HypoxiaBoundarySquared = 0;
	double QuiescenceBoundarySquared = 0;
	double TumourBoundarySquared = 0;
	double InnerBoundarySquared = 10000;
	double hypoxiaLevel = GetHypoxiaLevel();
	double quiescenceLevel = GetQuiescenceLevel();
	unsigned apoptoticCellCount = 0;
	unsigned tumourCellCount = 0;
	unsigned hypoxicTumourCellCount = 0;
	unsigned quiescentTumourCellCount = 0;
	unsigned normoxicTumourCellCount = 0;

	// Iterate over nodes, checking the cell type of each cell. If cell is a macrophage, write the location out
	for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
			node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
			++node_iter)
	{
		if (!node_iter->IsDeleted())
		{
			unsigned nodeIndex = node_iter->GetIndex();
			CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(nodeIndex);
			boost::shared_ptr<AbstractCellProliferativeType> celltype = pCell->GetCellProliferativeType();

			boost::shared_ptr<CellData> pCellData = pCell->GetCellData();
			double o2conc = pCellData->GetItem("oxygen");

			if (celltype->template IsType<MacrophageCellProliferativeType>())
			{

				const c_vector<double,DIM>& position = node_iter->rGetLocation();

				for (unsigned i=0; i<DIM; i++)
				{
					*pOutStreamMacrophageNodes << position[i] << " ";
				}
				*pOutStreamMacrophageOxygenConc << o2conc << " ";
			}

			if (pCell->HasApoptosisBegun())
			{
				apoptoticCellCount++;
			}

			else if (celltype->template IsType<StemCellProliferativeType>())
			{
				tumourCellCount++;
				const c_vector<double,DIM>& tumourCellLocation = node_iter->rGetLocation();

				// Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
				c_vector<double,DIM> displacement;
				displacement = centroid - tumourCellLocation;

				double distanceSquared = 0;
				for (unsigned i=0; i<DIM; i++)
				{
					distanceSquared += displacement[i]*displacement[i];
				}
				if (InnerBoundarySquared > distanceSquared)
				{
					InnerBoundarySquared = distanceSquared;
				}

				if (o2conc < hypoxiaLevel)
				{
					hypoxicTumourCellCount++;
					if (HypoxiaBoundarySquared < distanceSquared)
					{
						HypoxiaBoundarySquared = distanceSquared;
					}
				}
				else if (o2conc < quiescenceLevel && o2conc > hypoxiaLevel)
				{
					quiescentTumourCellCount++;
					if (QuiescenceBoundarySquared < distanceSquared)
					{
						QuiescenceBoundarySquared = distanceSquared;
					}
				}
				else if (o2conc > quiescenceLevel)
				{
					normoxicTumourCellCount++;
					if (TumourBoundarySquared < distanceSquared)
					{
						TumourBoundarySquared = distanceSquared;
					}
				}
			}

		}
	}
	HypoxiaBoundary = sqrt(HypoxiaBoundarySquared);
	QuiescenceBoundary = sqrt(QuiescenceBoundarySquared);
	TumourBoundary = sqrt(TumourBoundarySquared);
	InnerBoundary = sqrt(InnerBoundarySquared);
	*pOutStreamCompartmentBoundaries << InnerBoundary  << "\t" << HypoxiaBoundary << "\t" << QuiescenceBoundary << "\t" << TumourBoundary << "\t";

	*pOutStreamMacrophageNodes << std::endl;
	*pOutStreamCompartmentBoundaries << std::endl;
	*pOutStreamMacrophageOxygenConc << std::endl;
	*pOutStreamTumourCellCount << tumourCellCount << "\t" << apoptoticCellCount << "\t" << hypoxicTumourCellCount << "\t" << quiescentTumourCellCount << "\t" << normoxicTumourCellCount << std::endl;

}

template<unsigned DIM>
void OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Make sure the cell population is updated
	rCellPopulation.Update();

	/**
	 * We want to output these things:
	 *
	 * 1) All Macrophage node locations - This enables us to examine the distribution more closely during post-processing
	 * 2) The "size" of each tumour compartment: where is the centroid? Where do the hypoxic and quiescent sections begin?
	 * 3) Oxygen concentration at each macrophage
	 *
	 */
	double timestepDuration = SimulationTime::Instance()->GetTimeStep();


	std::string outputDirectory = GetOutputDirectory();

	double outputFrequency = GetOutputFrequencyInHours();
	//double time = SimulationTime::Instance()->GetTime();
	bool shouldWePrint = (timestepDuration*mTimestepsElapsed >= outputFrequency);
	if(shouldWePrint)
	//if(std::fmod(time,outputFrequency)==0)
	{
		//OutputFileHandler output_file_handler(outputDirectory+"/SummaryStatistics/", false);
		OutputFileHandler output_file_handler(outputDirectory, false);
		out_stream pOutStreamMacrophageNodes  = output_file_handler.OpenOutputFile("results.vizMacrophageNodes",std::ios::app);
		out_stream pOutStreamCompartmentBoundaries  = output_file_handler.OpenOutputFile("results.tumourCompartmentRadii",std::ios::app);
		out_stream pOutStreamMacrophageOxygenConc  = output_file_handler.OpenOutputFile("results.macrophageO2Concentration",std::ios::app);
		out_stream pOutStreamTumourCellCount  = output_file_handler.OpenOutputFile("results.tumourCellCount",std::ios::app);

		// We manually write the data to each of the files, rather than use the standard CellWriter scheme, as we want to bypass it so that no .vtu files are output in the NodeBased simulation
		WriteData(rCellPopulation, pOutStreamMacrophageNodes, pOutStreamCompartmentBoundaries, pOutStreamMacrophageOxygenConc, pOutStreamTumourCellCount);
	}
	mTimestepsElapsed++;

	if(shouldWePrint)
	{
		mTimestepsElapsed=1;
	}




}

template<unsigned DIM>
void OutputOnlyMacrophageSummaryStatisticsModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
	// No parameters to output, so just call method on direct parent class
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class OutputOnlyMacrophageSummaryStatisticsModifier<1>;
template class OutputOnlyMacrophageSummaryStatisticsModifier<2>;
template class OutputOnlyMacrophageSummaryStatisticsModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OutputOnlyMacrophageSummaryStatisticsModifier)


///*
//
//Copyright (c) 2005-2017, University of Oxford.
//All rights reserved.
//
//University of Oxford means the Chancellor, Masters and Scholars of the
//University of Oxford, having an administrative office at Wellington
//Square, Oxford OX1 2JD, UK.
//
//This file is part of Chaste.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of the University of Oxford nor the names of its
//   contributors may be used to endorse or promote products derived from this
//   software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//*/
//
//#include "NodeBasedCellPopulationWithVaryingDamping.hpp"
//
//
//template<unsigned DIM>
//NodeBasedCellPopulationWithVaryingDamping<DIM>::NodeBasedCellPopulationWithVaryingDamping(NodesOnlyMesh<DIM>& rMesh,
//                                      std::vector<CellPtr>& rCells,
//                                      const std::vector<unsigned> locationIndices,
//                                      bool deleteMesh)
//    : NodeBasedCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh, false),
//	  mBoundaryNodeDampingConstant(1.0),
//	  mInteriorNodeDampingConstant(1.0)
//{
//}
//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//double AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetBoundaryNodeDampingConstant()
//{
//	return mBoundaryNodeDampingConstant;
//}
//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//void AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetBoundaryNodeDampingConstant(double newBoundaryNodeDampingConstant)
//{
//	mBoundaryNodeDampingConstant = newBoundaryNodeDampingConstant;
//}
//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//double AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetInteriorNodeDampingConstant()
//{
//	return mInteriorNodeDampingConstant;
//}
//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//void AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetInteriorNodeDampingConstant(double newInteriorNodeDampingConstant)
//{
//	mInteriorNodeDampingConstant = newInteriorNodeDampingConstant;
//}
//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//double AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetDampingConstant(unsigned nodeIndex)
//{
//	// Overwritten GetDampingConstant method
//	// If our node is on the boundary, we return mBoundaryNodeDampingConstant. Otherwise, we return mInteriorNodeDampingConstant
//
//	if(this->GetNode(nodeIndex)->IsBoundaryNode())
//	{
//		return mBoundaryNodeDampingConstant;
//	}
//	else
//	{
//		return mInteriorNodeDampingConstant;
//	}
//
//}
//
//template<unsigned DIM>
//void NodeBasedCellPopulationWithVaryingDamping<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
//{
//	*rParamsFile << "\t\t<mBoundaryNodeDampingConstant>" << mBoundaryNodeDampingConstant << "</mBoundaryNodeDampingConstant>\n";
//	*rParamsFile << "\t\t<mInteriorNodeDampingConstant>" << mInteriorNodeDampingConstant << "</mInteriorNodeDampingConstant>\n";
//
//    // Call method on direct parent class
//    NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
//}
//
//// Explicit instantiation
//template class NodeBasedCellPopulationWithVaryingDamping<1>;
//template class NodeBasedCellPopulationWithVaryingDamping<2>;
//template class NodeBasedCellPopulationWithVaryingDamping<3>;
//
//// Serialization for Boost >= 1.36
//#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithVaryingDamping)

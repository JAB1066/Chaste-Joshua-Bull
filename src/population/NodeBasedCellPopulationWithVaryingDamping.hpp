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
//#ifndef NODEBASEDCELLPOPULATIONWITHVARYINGDAMPING_HPP_
//#define NODEBASEDCELLPOPULATIONWITHVARYINGDAMPING_HPP_
//
//#include "NodeBasedCellPopulation.hpp"
//#include "NodesOnlyMesh.hpp"
//
//#include "ChasteSerialization.hpp"
//#include <boost/serialization/base_object.hpp>
//
///**
// * A NodeBasedCellPopulationWithVaryingDamping is a NodeBasedCellPopulation with the ability to set different damping coefficients inside and outside the tumour.
// * This permits simulation of the spheroid growing in to a medium which provides more resistance than the area already occupied by the tumour
// */
//template<unsigned DIM>
//class NodeBasedCellPopulationWithVaryingDamping : public NodeBasedCellPopulation<DIM>
//{
//
//protected:
//    /**
//     * Check consistency of our internal data structures.
//     */
//    void Validate();
//
//    /**
//     * Overridden AcceptCellWritersAcrossPopulation() method.
//     *
//     * Calls #AcceptCellWriter() across the whole population,
//     * iterating in an appropriate way to skip particle nodes.
//     */
//    virtual void AcceptCellWritersAcrossPopulation();
//
//private:
//
//    double mBoundaryNodeDampingConstant;
//
//    double mInteriorNodeDampingConstant;
//
//    /** Needed for serialization. */
//    friend class boost::serialization::access;
//    /**
//     * Serialize the object and its member variables.
//     *
//     * Note that serialization of the nodes is handled by load/save_construct_data,
//     * so we don't actually have to do anything here except delegate to the base class.
//     *
//     * @param archive the archive
//     * @param version the current version of this class
//     */
//    template<class Archive>
//    void serialize(Archive & archive, const unsigned int version)
//    {
//        archive & boost::serialization::base_object<NodeBasedCellPopulation<DIM> >(*this);
//    }
//
//public:
//
//    /**
//     * Default constructor.
//     *
//     * Note that the cell population will take responsibility for freeing the memory used by the nodes.
//     *
//     * @param rMesh a mutable nodes-only mesh
//     * @param rCells a vector of cells
//     * @param locationIndices an optional vector of location indices that correspond to real cells
//     * @param deleteMesh whether to delete nodes-only mesh in destructor
//     */
//    NodeBasedCellPopulationWithVaryingDamping(NodesOnlyMesh<DIM>& rMesh,
//                            std::vector<CellPtr>& rCells,
//                            const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
//                            bool deleteMesh=false);
//
//
//    double GetBoundaryNodeDampingConstant();
//
//    void SetBoundaryNodeDampingConstant(double newBoundaryNodeDampingConstant);
//
//    double GetInteriorNodeDampingConstant();
//
//    void SetInteriorNodeDampingConstant(double newInteriorNodeDampingConstant);
//
//
//    /**
//     * Outputs CellPopulation parameters to file
//     *
//     * As this method is pure virtual, it must be overridden
//     * in subclasses.
//     *
//     * @param rParamsFile the file stream to which the parameters are output
//     */
//    void OutputCellPopulationParameters(out_stream& rParamsFile);
//};
//
//#include "SerializationExportWrapper.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithVaryingDamping)
//
//namespace boost
//{
//namespace serialization
//{
///**
// * Serialize information required to construct a NodeBasedCellPopulationWithVaryingDamping.
// */
//template<class Archive, unsigned DIM>
//inline void save_construct_data(
//    Archive & ar, const NodeBasedCellPopulationWithVaryingDamping<DIM> * t, const unsigned int file_version)
//{
//    // Save data required to construct instance
//    const NodesOnlyMesh<DIM>* p_mesh = &(t->rGetMesh());
//    ar & p_mesh;
//}
//
///**
// * De-serialize constructor parameters and initialise a NodeBasedCellPopulationWithParticles.
// * Loads the mesh from separate files.
// */
//template<class Archive, unsigned DIM>
//inline void load_construct_data(
//    Archive & ar, NodeBasedCellPopulationWithVaryingDamping<DIM> * t, const unsigned int file_version)
//{
//    // Retrieve data from archive required to construct new instance
//    NodesOnlyMesh<DIM>* p_mesh;
//    ar >> p_mesh;
//
//    // Invoke inplace constructor to initialise instance
//    ::new(t)NodeBasedCellPopulationWithVaryingDamping<DIM>(*p_mesh);
//}
//}
//} // namespace ...
//
//#endif /*NODEBASEDCELLPOPULATIONWITHVARYINGDAMPING_HPP_*/

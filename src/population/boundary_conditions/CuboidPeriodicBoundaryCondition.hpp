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

#ifndef CUBOIDPERIODICBOUNDARYCONDITION_HPP_
#define CUBOIDPERIODICBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "ChastePoint.hpp"

#include "Debug.hpp"
/**
 * A boundary condition class for a node based cell population within a cuboid. When
 * the force on a node causes it to pass through a face of the cuboid, it reenters
 * through the opposite face at the corresponding location.
 *
 * Note: Upper and Lower corners of the cuboid are assumed to have upper[i] > lower[i] for all i
 */
template<unsigned DIM>
class CuboidPeriodicBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /** The lower corner of the cuboid */
    ChastePoint<DIM> mLower;

    /** The upper corner of the cuboid */
    ChastePoint<DIM> mUpper;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
        archive & mLower;
        archive & mUpper;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param lower the lower corner of the cuboid
     * @param upper the upper corner of the cuboid
     */
    CuboidPeriodicBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                    ChastePoint<DIM> lower,
									ChastePoint<DIM> upper);

    /**
     * @return #mUpper.
     */
    ChastePoint<DIM> GetLowerPoint() const;

    /**
     * @return #mLower.
     */
    ChastePoint<DIM> GetUpperPoint() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};
//
//#include "SerializationExportWrapper.hpp"
//CHASTE_CLASS_EXPORT(CuboidPeriodicBoundaryCondition)
//#include "SerializationExportWrapperForCpp.hpp"
//CHASTE_CLASS_EXPORT(CuboidPeriodicBoundaryCondition)


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CuboidPeriodicBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CuboidPeriodicBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CuboidPeriodicBoundaryCondition<DIM>* t, const unsigned int file_version)
{

    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive chaste points one component at a time
    c_vector<double, DIM> lower = t->GetLowerPoint().rGetLocation();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << lower[i];
    }
    c_vector<double, DIM> upper = t->GetUpperPoint().rGetLocation();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << upper[i];
    }

}

/**
 * De-serialize constructor parameters and initialize a CuboidPeriodicBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CuboidPeriodicBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, DIM> lower_vec;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> lower_vec[i];
    }
    c_vector<double, DIM> upper_vec;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> upper_vec[i];
    }

    ChastePoint<DIM> lower (lower_vec);
    ChastePoint<DIM> upper (upper_vec);
    // Invoke inplace constructor to initialise instance
    ::new(t)CuboidPeriodicBoundaryCondition<DIM>(p_cell_population, lower, upper);
}
}
} // namespace ...


#endif /*CUBOIDPERIODICBOUNDARYCONDITION_HPP_*/

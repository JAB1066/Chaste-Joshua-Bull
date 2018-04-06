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

#ifndef UNIFORMCELLCYCLEMODELWITHQUIESCENCE_HPP_
#define UNIFORMCELLCYCLEMODELWITHQUIESCENCE_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A stochastic cell-cycle model where cells divide with a stochastic cell cycle duration
 * with the length of the cell cycle drawn from a uniform distribution
 * on [mMinCellCycleDuration, mMaxCellCycleDuration].
 *
 * If the cell is differentiated, then the cell cycle duration is set to be infinite,
 * so that the cell will never divide.
 */
class UniformCellCycleModelWithQuiescence : public AbstractSimpleCellCycleModel
{
	friend class TestSimpleCellCycleModels;

private:

	/**
	 * The minimum cell cycle duration. Used to define the uniform distribution.
	 * Defaults to 12 hours.
	 */
	double mMinCellCycleDuration;

	/**
	 * The maximum cell cycle duration. Used to define the uniform distribution.
	 * Defaults to 14 hours.
	 */
	double mMaxCellCycleDuration;

	/** Needed for serialization. */
	friend class boost::serialization::access;
	/**
	 * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
	 *
	 * @param archive the archive
	 * @param version the current version of this class
	 */
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
		archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

		// Make sure the RandomNumberGenerator singleton gets saved too
		SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
		archive & p_wrapper;
		archive & mMinCellCycleDuration;
		archive & mMaxCellCycleDuration;
		archive & mCurrentHypoxicDuration;
		archive & mHypoxicConcentration;
		archive & mQuiescentConcentration;
		archive & mCriticalHypoxicDuration;
	}

protected:

	/**
	 * How long the current period of hypoxia has lasted.
	 * Has units of hours.
	 */
	double mCurrentHypoxicDuration;

	/**
	 * Non-dimensionalized oxygen concentration below which cells are
	 * considered to be hypoxic. A prolonged period of hypoxia causes
	 * the cell to become apoptotic.
	 */
	double mHypoxicConcentration;

	/**
	 * Non-dimensionalized oxygen concentration below which cells are
	 * considered to be quiescent and slow their progress through the
	 * G1 phase of the cell cycle.
	 */
	double mQuiescentConcentration;

	/**
	 * Non-dimensionalized critical hypoxic duration.
	 * Has units of hours.
	 */
	double mCriticalHypoxicDuration;

    /**
     * The time when the current period of hypoxia began.
     */
    double mCurrentHypoxiaOnsetTime;

	/**
	 * Protected copy-constructor for use by CreateCellCycleModel().
	 *
	 * The only way for external code to create a copy of a cell cycle model
	 * is by calling that method, to ensure that a model of the correct subclass is created.
	 * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
	 *
	 * This method is called by child classes to set member variables for a daughter cell upon cell division.
	 * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
	 * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
	 * can be done in InitialiseDaughterCell().
	 *
	 * @param rModel the cell cycle model to copy.
	 */
    UniformCellCycleModelWithQuiescence(const UniformCellCycleModelWithQuiescence& rModel);

public:

	/**
	 * Constructor - just a default, mBirthTime is set in the AbstractCellCycleModel class.
	 */
    UniformCellCycleModelWithQuiescence();

	/**
	 * Overridden SetCellCycleDuration() method to add stochastic cell cycle times
	 */
	void SetCellCycleDuration();

	/**
	 * Overridden builder method to create new copies of
	 * this cell-cycle model.
	 *
	 * @return new cell-cycle model
	 */
	AbstractCellCycleModel* CreateCellCycleModel();

	/**
	 * @return mMinCellCycleDuration
	 */
	double GetMinCellCycleDuration();

	/**
	 * Set mMinCellCycleDuration.
	 *
	 * @param minCellCycleDuration
	 */
	void SetMinCellCycleDuration(double minCellCycleDuration);

	/**
	 * @return mMaxCellCycleDuration
	 */
	double GetMaxCellCycleDuration();

	/**
	 * Set mMaxCellCycleDuration.
	 *
	 * @param maxCellCycleDuration
	 */
	void SetMaxCellCycleDuration(double maxCellCycleDuration);

	/**
	 * Overridden GetAverageTransitCellCycleTime() method.
	 *
	 * @return the average of mMinCellCycleDuration and mMaxCellCycleDuration
	 */
	double GetAverageTransitCellCycleTime();

	/**
	 * Overridden GetAverageStemCellCycleTime() method.
	 *
	 * @return the average of mMinCellCycleDuration and mMaxCellCycleDuration
	 */
	double GetAverageStemCellCycleTime();

    /**
     * @return mCurrentHypoxicDuration
     */
    double GetCurrentHypoxicDuration() const;

    /**
      * @return mHypoxicConcentration
      */
     double GetHypoxicConcentration() const;

     /**
      * Set method for mHypoxicConcentration.
      *
      * @param hypoxicConcentration the new value of mHypoxicConcentration
      */
     void SetHypoxicConcentration(double hypoxicConcentration);

     /**
      * @return mQuiescentConcentration
      */
     double GetQuiescentConcentration() const;

     /**
      * Set method for mQuiescentConcentration.
      *
      * @param quiescentConcentration the new value of mQuiescentConcentration
      */
     void SetQuiescentConcentration(double quiescentConcentration);

     /**
      * @return mCriticalHypoxicDuration
      */
     double GetCriticalHypoxicDuration() const;

     /**
      * Set method for mCriticalHypoxicDuration.
      *
      * @param criticalHypoxicDuration the new value of mCriticalHypoxicDuration
      */
     void SetCriticalHypoxicDuration(double criticalHypoxicDuration);

     /**
      * Set method for mCurrentHypoxiaOnsetTime.
      *
      * @param currentHypoxiaOnsetTime the new value of mCurrentHypoxiaOnsetTime
      */
     void SetCurrentHypoxiaOnsetTime(double currentHypoxiaOnsetTime);

     void UpdateHypoxicDuration();

     bool ReadyToDivide();

	/**
	 * Overridden OutputCellCycleModelParameters() method.
	 *
	 * @param rParamsFile the file stream to which the parameters are output
	 */
	virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(UniformCellCycleModelWithQuiescence)

#endif /*UNIFORMCELLCYCLEMODELWITHQUIESCENCE_HPP_*/

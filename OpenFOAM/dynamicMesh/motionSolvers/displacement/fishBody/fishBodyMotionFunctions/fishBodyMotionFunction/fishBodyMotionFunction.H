/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::fishBodyMotionFunctions

Description
    Namespace for solid-body motions


Class
    Foam::fishBodyMotionFunction

Description
    Base class for defining solid-body motions

SourceFiles
    fishBodyMotionFunction.C
    dynamicFvMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef fishBodyMotionFunction_H
#define fishBodyMotionFunction_H

#include "Time.H"
#include "dictionary.H"
#include "septernion.H"
#include "pointFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fishBodyMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class fishBodyMotionFunction
{
protected:

    // Protected Data

        dictionary SBMFCoeffs_;

        const Time& time_;

        pointField p0_;

    // Protected Member Functions

        //- No copy construct
        fishBodyMotionFunction(const fishBodyMotionFunction&) = delete;

        //- No copy assignment
        void operator=(const fishBodyMotionFunction&) = delete;


public:

    //- Runtime type information
    TypeName("fishBodyMotionFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fishBodyMotionFunction,
            dictionary,
            (const dictionary& SBMFCoeffs, const Time& runTime),
            (SBMFCoeffs, runTime)
        );


    // Constructors

        //- Construct from the SBMFCoeffs dictionary and Time
        fishBodyMotionFunction
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<fishBodyMotionFunction> clone() const = 0;


    // Selectors

        //- Select constructed from the SBMFCoeffs dictionary and Time
        static autoPtr<fishBodyMotionFunction> New
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~fishBodyMotionFunction() = default;


    // Member Functions

        //- Return the transformed pointField
        virtual tmp<vectorField> transformationPoints(pointField& p0_) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) = 0;

        //- Write in dictionary format
        virtual void writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

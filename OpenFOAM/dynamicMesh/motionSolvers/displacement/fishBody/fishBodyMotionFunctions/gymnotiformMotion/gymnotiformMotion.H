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

Class
    Foam::fishBodyMotionFunctions::gymnotiformMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating rotation.

SourceFiles
    gymnotiformMotion.C

\*---------------------------------------------------------------------------*/

#ifndef gymnotiformMotion_H
#define gymnotiformMotion_H

#include "fishBodyMotionFunction.H"
#include "primitiveFields.H"
#include "pointFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fishBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class gymnotiformMotion Declaration
\*---------------------------------------------------------------------------*/

class gymnotiformMotion
:
    public fishBodyMotionFunction
{
    // Private Data

        //- Centre of gravity
        point origin_;

        //- Amplitude
        scalar amplitude_;

        //- Wave number
        scalar waveNumber_;

        //- Length
        scalar length_;

        //- Ramp
        scalar ramp_;

        //- Radial velocity
        scalar omega_;

        //- Max Angle
        scalar maxAngle_;

        //- Define a pointField
        pointField p0_;

    // Private Member Functions

        //- No copy construct
        gymnotiformMotion(const gymnotiformMotion&) = delete;

        //- No copy assignment
        void operator=(const gymnotiformMotion&) = delete;


public:

    //- Runtime type information
    TypeName("gymnotiformMotion");


    // Constructors

        //- Construct from components
        gymnotiformMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<fishBodyMotionFunction> clone() const
        {
            return autoPtr<fishBodyMotionFunction>
            (
                new gymnotiformMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~gymnotiformMotion() = default;


    // Member Functions

        //- Return the transformed pointField
        virtual tmp<vectorField> transformationPoints(pointField& p0_) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fishBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

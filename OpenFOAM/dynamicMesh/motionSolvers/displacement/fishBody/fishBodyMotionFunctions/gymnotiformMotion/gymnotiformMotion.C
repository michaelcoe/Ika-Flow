/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "gymnotiformMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fishBodyMotionFunctions
{
    defineTypeNameAndDebug(gymnotiformMotion, 0);
    addToRunTimeSelectionTable
    (
        fishBodyMotionFunction,
        gymnotiformMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fishBodyMotionFunctions::gymnotiformMotion::
gymnotiformMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    fishBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::pointField>
Foam::fishBodyMotionFunctions::gymnotiformMotion::
transformationPoints(pointField& p0) const
{
    scalar tm = time_.value();

    forAll(p0, pointI)
    {
    	const scalar x = (p0[pointI].component(0)-origin_[0])/length_;
        const scalar y = p0[pointI].component(1)-origin_[1];
        const scalar z = p0[pointI].component(2)-origin_[2];
                   
        const scalar yr = y + amplitude_ * maxAngle_ * sin(waveNumber_ * x - omega_ * tm) * length_;
        
        p0[pointI] = vector(x, yr, z);
    }

    return p0;
}


bool Foam::fishBodyMotionFunctions::gymnotiformMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    fishBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.readEntry("origin", origin_);
    SBMFCoeffs_.readEntry("amplitude", amplitude_);
    SBMFCoeffs_.readEntry("waveNumber", waveNumber_);
    SBMFCoeffs_.readEntry("length", length_);
    SBMFCoeffs_.readEntry("ramp", ramp_);
    SBMFCoeffs_.readEntry("omega", omega_);
    SBMFCoeffs_.readEntry("maxAngle", maxAngle_);

    return true;
}


// ************************************************************************* //

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

#include "ostraciiformMotion.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fishBodyMotionFunctions
{
    defineTypeNameAndDebug(ostraciiformMotion, 0);
    addToRunTimeSelectionTable
    (
        fishBodyMotionFunction,
        ostraciiformMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fishBodyMotionFunctions::ostraciiformMotion::
ostraciiformMotion
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
Foam::fishBodyMotionFunctions::ostraciiformMotion::
transformationPoints(pointField& p0) const
{
    const scalar tm = time_.value();

    forAll(p0, pointI)
    {
    	const scalar x = (p0[pointI].component(0)-origin_[0])/length_;
        const scalar y = p0[pointI].component(1)-origin_[1];
        const scalar z = p0[pointI].component(2)-origin_[2];

        scalar yr = 0;

        if (x >= pivot_)
        {
            // new value by equation
            const scalar xPivot = x - pivot_;
            
            const scalar thetaT = maxAngle_ * sin(waveNumber_*pivot_ - omega_*tm + phase_);

            yr = y + xPivot * tan(thetaT) * length_;

        }
        else
        {
            yr = y;
        }

        p0[pointI] = vector(x, yr, z);
    }

    return p0;
}


bool Foam::fishBodyMotionFunctions::ostraciiformMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    fishBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.readEntry("origin", origin_);
    SBMFCoeffs_.readEntry("waveNumber", waveNumber_);
    SBMFCoeffs_.readEntry("length", length_);
    SBMFCoeffs_.readEntry("ramp", ramp_);
    SBMFCoeffs_.readEntry("omega", omega_);
    SBMFCoeffs_.readEntry("pivot", pivot_);
    SBMFCoeffs_.readEntry("maxAngle", maxAngle_);
    SBMFCoeffs_.readEntry("phaseAngle", phase_);

    return true;
}


// ************************************************************************* //

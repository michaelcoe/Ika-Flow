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

\*---------------------------------------------------------------------------*/

#include "fishWavePointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fishMotionPointPatchVectorField::
fishMotionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    amplitude_(Zero),
    origin_(Zero),
    axis_(Zero),
    omega_(0.0),
    length_(0.0),
    waveNumber_(Zero)
{}


Foam::fishMotionPointPatchVectorField::
fishMotionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    amplitude_(dict.lookup("amplitude")),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    omega_(dict.get<scalar>("omega")),
    length_(dict.get<scalar>("length")),
    waveNumber_(dict.getOrDefault<vector>("waveNumber", Zero))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::fishMotionPointPatchVectorField::
fishMotionPointPatchVectorField
(
    const fishMotionPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    omega_(ptf.omega_),
    length_(ptf.length_),
    waveNumber_(ptf.waveNumber_)
{}


Foam::fishMotionPointPatchVectorField::
fishMotionPointPatchVectorField
(
    const fishMotionPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    amplitude_(ptf.amplitude_),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    omega_(ptf.omega_),
    length_(ptf.length_),
    waveNumber_(ptf.waveNumber_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fishMotionPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    const vectorField p0Rel(patch().localPoints() - origin_);
    const scalarField points( waveNumber_ & p0Rel);
    const scalarField xCoord = p0Rel.component(vector::X)/length_;

    Field<vector>::operator=
    (
        axis_ *
        (amplitude_[0] + (amplitude_[1]*xCoord) + (amplitude_[2]*xCoord*xCoord))
        *sin(omega_*t.value() - points)*length_
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::fishMotionPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("amplitude", amplitude_);
    os.writeEntry("origin", origin_);
    os.writeEntry("axis", axis_);
    os.writeEntry("omega", omega_);
    os.writeEntry("length", length_);
    os.writeEntry("waveNumber", waveNumber_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        fishMotionPointPatchVectorField
    );
}

// ************************************************************************* //

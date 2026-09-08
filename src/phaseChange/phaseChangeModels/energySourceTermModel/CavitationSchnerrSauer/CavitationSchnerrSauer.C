/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "CavitationSchnerrSauer.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CavitationSchnerrSauer, 0);
    addToRunTimeSelectionTable
    (
        energySourceTermModel,
        CavitationSchnerrSauer,
        components
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CavitationSchnerrSauer::CavitationSchnerrSauer
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const compressibleInterPhaseTransportModel& turbModel,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    reconstructionSchemes& surf,
    const dictionary& dict
)
:
    energySourceTermModel
    (
        typeName,
        phase1,
        phase2,
        turbModel,
        p,
        satModel,
        surf,
        dict
    )
{

}

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix> Foam::CavitationSchnerrSauer::TSource1()
{
    // Cavitation model does not directly contribute to phase1 temperature equation
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            phase1_.thermo().T(),
            dimEnergy/dimTime
        )
    );
}

Foam::tmp<Foam::fvScalarMatrix> Foam::CavitationSchnerrSauer::TSource2()
{
    // Cavitation model does not directly contribute to phase2 temperature equation
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            phase2_.thermo().T(),
            dimEnergy/dimTime
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::CavitationSchnerrSauer::energySource()
{
    return volScalarField::New
    (
        "zeroEnergySource",
        phase1_.mesh(),
        dimensionedScalar(dimPower/dimVolume, 0)
    );
}

Foam::tmp<Foam::volScalarField> Foam::CavitationSchnerrSauer::energyFlux1()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "energyFlux1",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phase1_.mesh(),
            dimensionedScalar(dimPower/dimArea, 0)
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::CavitationSchnerrSauer::energyFlux2()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "energyFlux2",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phase1_.mesh(),
            dimensionedScalar(dimPower/dimArea, 0)
        )
    );
}

// ************************************************************************* //

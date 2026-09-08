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

#include "HybridPhaseChangeMass.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HybridPhaseChangeMass, 0);
    addToRunTimeSelectionTable
    (
        massSourceTermModel,
        HybridPhaseChangeMass,
        components
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HybridPhaseChangeMass::HybridPhaseChangeMass
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    reconstructionSchemes& surf,
    const dictionary& dict
)
:
    massSourceTermModel
    (
        typeName,
        phase1,
        phase2,
        p,
        satModel,
        surf,
        dict
    ),
    interfaceModel_(),
    cavitationModel_()
{
    // 修复：不需要额外传入 "type" 字符串
    // New 函数会从字典中读取 "type" 条目
    const dictionary& interfaceDict = modelDict().subDict("interfaceModel");
    interfaceModel_ = massSourceTermModel::New
    (
        phase1_,
        phase2_,
        p_,
        satModel_,
        surf_,
        interfaceDict
    );
    
    const dictionary& cavitationDict = modelDict().subDict("cavitationModel");
    cavitationModel_ = massSourceTermModel::New
    (
        phase1_,
        phase2_,
        p_,
        satModel_,
        surf_,
        cavitationDict
    );
}

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::HybridPhaseChangeMass::zeroInterfaceCells
(
    volScalarField& field
) const
{
    const volScalarField& alpha1 = phase1_;
    const scalar cutOff = 1e-3;

    forAll(alpha1, cellI)
    {
        if
        (
            alpha1[cellI] > cutOff
         && alpha1[cellI] < (1.0 - cutOff)
        )
        {
            field[cellI] = 0.0;
        }
    }
}

// * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> 
Foam::HybridPhaseChangeMass::massSource(volScalarField& rhoSource)
{
    return interfaceModel_->massSource(rhoSource);
}

Foam::tmp<Foam::volScalarField>
Foam::HybridPhaseChangeMass::alphaSource( volScalarField& rhoSource)
{
    tmp<volScalarField> alphaSource
    (
        rhoSource / phase1_.thermo().rho()
    );

   return alphaSource;
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::HybridPhaseChangeMass::vDotAlphal() const
{
    Pair<tmp<volScalarField>> vDot = cavitationModel_->vDotAlphal();

    volScalarField& vDotc = vDot[0].ref();
    volScalarField& vDotv = vDot[1].ref();

    zeroInterfaceCells(vDotc);
    zeroInterfaceCells(vDotv);

    return vDot;
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::HybridPhaseChangeMass::vDotP() const
{
    Pair<tmp<volScalarField>> vDot = cavitationModel_->vDotP();

    volScalarField& vDotc = vDot[0].ref();
    volScalarField& vDotv = vDot[1].ref();

    zeroInterfaceCells(vDotc);
    zeroInterfaceCells(vDotv);

    return vDot;
}
// ************************************************************************* //

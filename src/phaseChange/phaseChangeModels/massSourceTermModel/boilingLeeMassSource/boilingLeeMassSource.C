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
-------------------------------------------------------------------------------
            Copyright (c) 2017-2024, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
            Copyright (c) 2025-2028, University College Dublin (UCD)
-------------------------------------------------------------------------------
Class
    Foam::boilingLeeMassSource

Description
    Volumetric Lee evaporation/condensation model.
    BoilingLee is the energySourceTerm, while this scripts, 
    boilingLeeMassSource is the massSourceTerm associated to the same Lee model.

SourceFiles
    boilingLeeMassSource.C
\*---------------------------------------------------------------------------*/

#include "boilingLeeMassSource.H"
#include "addToRunTimeSelectionTable.H"
//#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
//#include "fvcSnGrad.H"

//#include "centredCPCCellToCellStencilObject.H"
//#include "SortableList.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(boilingLeeMassSource, 0);
    addToRunTimeSelectionTable(massSourceTermModel,boilingLeeMassSource,components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::boilingLeeMassSource::boilingLeeMassSource
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
    C_(modelDict().get<dimensionedScalar>("C")),
    Tactivate_(modelDict().get<dimensionedScalar>("Tactivate")),
    alphaMin_(modelDict().lookupOrDefault<scalar>("alphaMin",0.0))
{}

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //
// ------------------------------------------------------------------------- //
// Full explicit Lee mass-transfer rate
// phase1 = liquid
// mDot = C * alpha1 * rho1 * (T - Tactivate)/Tactivate 
// @ C > 0: evaporation when T > Tactivate
// Dimensions: kg/(m3 s)
// ------------------------------------------------------------------------- //

Foam::tmp<Foam::volScalarField>Foam::boilingLeeMassSource::Kexp(const volScalarField& T) const
{
    const volScalarField& rho1 = phase1_.thermo().rho();

    volScalarField alpha1Limited
    (
        min(max(phase1_, scalar(0)), scalar(1))
    );

    volScalarField coeff
    (
        C_
       *alpha1Limited
       *rho1
       *pos(alpha1Limited - alphaMin_)
       *(T - Tactivate_)
       /Tactivate_
    );

    if (sign(C_.value()) > 0)
    {
        return coeff*pos(T - Tactivate_);
    }
    else
    {
        return coeff*pos(Tactivate_ - T);
    }
}

// * * * * * * * * * * Public Member Functions * * * * * * * * * //

//Foam::tmp<Foam::volScalarField>Foam::boilingLeeMassSource::alphaTransfer(const volScalarField& T) const
//{
//    const volScalarField& rho1 = phase1_.thermo().rho();
//    const volScalarField& rho2 = phase2_.thermo().rho();
//    const volScalarField& alpha1 = phase1_;
//
//    tmp<volScalarField> tmDot = Kexp(T);
//    const volScalarField& mDot = tmDot();
//
//    return
//       -mDot
//       *(
//            1.0/rho1
//          - alpha1*(1.0/rho1 - 1.0/rho2)
//        );
//}

Foam::tmp<Foam::volScalarField>Foam::boilingLeeMassSource::dmdt(const volScalarField& T) const
{
    // Sign convention consistent with icoReacting:
    //      dmdt > 0 : mass transfer INTO phase1
    //      dmdt < 0 : mass transfer OUT OF phase1
    //
    // phase1 = liquid
    // Kexp > 0 = evaporation
    //
    // Therefore evaporation gives dmdt < 0.

    return -Kexp(T);
}

Foam::tmp<Foam::volScalarField>Foam::boilingLeeMassSource::volTransfer(const volScalarField& T) const
{
    const volScalarField& rho1 = phase1_.thermo().rho();
    const volScalarField& rho2 = phase2_.thermo().rho();

    tmp<volScalarField> tdmdt = dmdt(T);

    return
        tdmdt()
       *(1.0/rho1 - 1.0/rho2);
}

Foam::tmp<Foam::volScalarField> Foam::boilingLeeMassSource::massSource(volScalarField& rhoSource)
{
    return rhoSource * 0.0;
}

Foam::tmp<Foam::volScalarField> Foam::boilingLeeMassSource::alphaSource(volScalarField& rhoSource)
{
    return rhoSource / phase1_.thermo().rho() * 0.0;
}

// ************************************************************************* //

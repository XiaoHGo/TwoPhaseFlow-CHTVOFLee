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
    Foam::BoilingLee

Description
    Volumetric Lee evaporation/condensation model.
    The model follows the source-term decomposition used by 
    MassTransferPhaseSystem:
        mDot = KSp(variable)*variable + KSu(variable)
    together with a full explicit source
        mDot = Kexp(T,p,...)
    This interface allows the source-term architecture to support 
    temperature-driven, pressure-driven, or mixed phase-change models.
    BoilingLee itself is presently temperature-driven.

SourceFiles
    BoilingLee.C
\*---------------------------------------------------------------------------*/

#include "BoilingLee.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BoilingLee, 0);
    addToRunTimeSelectionTable
    (
        energySourceTermModel,
        BoilingLee,
        components
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BoilingLee::BoilingLee
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
    ),
    C_(modelDict().get<dimensionedScalar>("C")),
    Tactivate_(modelDict().get<dimensionedScalar>("Tactivate")),
    alphaMin_(modelDict().lookupOrDefault<scalar>("alphaMin",0.0))
{

}

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
// Full explicit Lee mass-transfer rate
Foam::tmp<Foam::volScalarField>Foam::BoilingLee::Kexp(const volScalarField& T) const
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
        return coeff * pos(T - Tactivate_);
    }
    else
    {
        // Retained for consistency with generic Lee.
        // For real condensation the donor phase should normally be vapour.
        return coeff * pos(Tactivate_ - T);
    }
}


// Temperature-linearised implicit coefficient
Foam::tmp<Foam::volScalarField>Foam::BoilingLee::KSp(const volScalarField& T) const
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
       /Tactivate_
    );

    if (sign(C_.value()) > 0)
    {
        return coeff*pos(T - Tactivate_);
    }
    else
    {
        return -coeff*pos(Tactivate_ - T);
    }
}


// Explicit part of Lee linearisation
Foam::tmp<Foam::volScalarField>Foam::BoilingLee::KSu(const volScalarField& T) const
{
    const volScalarField& rho1 = phase1_.thermo().rho();

    volScalarField alpha1Limited(min(max(phase1_, scalar(0)), scalar(1)));

    volScalarField coeff
    (
        C_
       *alpha1Limited
       *rho1
       *pos(alpha1Limited - alphaMin_)
    );

    if (sign(C_.value()) > 0)
    {
        return -coeff*pos(T - Tactivate_);
    }
    else
    {
        return -coeff*pos(Tactivate_ - T);
    }
}


// Baseline: no additional spatial redistribution
//Foam::tmp<Foam::volScalarField>Foam::BoilingLee::smear(const volScalarField& source) const
//{
//    return scalar(1)*source;
//}


// Semi-implicit global-T latent heat source
Foam::tmp<Foam::fvScalarMatrix>Foam::BoilingLee::heatTransfer(const volScalarField& T)
{
    tmp<fvScalarMatrix> tEqn
    (
        new fvScalarMatrix
        (
            T,
            dimEnergy/dimTime
        )
    );

    fvScalarMatrix& eqn = tEqn.ref();

    tmp<volScalarField> tSp = KSp(T);
    tmp<volScalarField> tSu = KSu(T);

    const volScalarField& Sp = tSp();
    const volScalarField& Su = tSu();

    const volScalarField& L = satModel_.L();

    eqn -= fvm::Sp(Sp * L, T) + Su * L;

    return tEqn;
}


// Legacy two-T interfaces: zero for global-T Lee model
Foam::tmp<Foam::fvScalarMatrix>Foam::BoilingLee::TSource1()
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            phase1_.thermo().T(),
            dimEnergy/dimTime
        )
    );
}


Foam::tmp<Foam::fvScalarMatrix>Foam::BoilingLee::TSource2()
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            phase2_.thermo().T(),
            dimEnergy/dimTime
        )
    );
}


// Fully explicit source for diagnostics/post-processing
Foam::tmp<Foam::volScalarField>Foam::BoilingLee::energySource()
{
    //const volScalarField& T = phase1_.thermo().T();
    //tmp<volScalarField> tmDot = Kexp(T);
    //const volScalarField& L = satModel_.L();
    //tmp<volScalarField> tq
    //(
    //    -tmDot()*L
    //);
    //return smear(tq());
    //return tq();

     return volScalarField::New
    (
        "zeroEnergySource",
        phase1_.mesh(),
        dimensionedScalar(dimPower/dimVolume, 0)
    );
}


// No interface heat-flux treatment
Foam::tmp<Foam::volScalarField>Foam::BoilingLee::energyFlux1()
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
            dimensionedScalar(dimPower/dimArea,0)
        )
    );
}


Foam::tmp<Foam::volScalarField>Foam::BoilingLee::energyFlux2()
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
            dimensionedScalar(dimPower/dimArea,0)
        )
    );
}

// ************************************************************************* //

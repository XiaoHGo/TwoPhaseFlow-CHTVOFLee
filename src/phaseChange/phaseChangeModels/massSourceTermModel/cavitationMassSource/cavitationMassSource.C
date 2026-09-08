/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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


#include "cavitationMassSource.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

#include "centredCPCCellToCellStencilObject.H"
#include "SortableList.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cavitationMassSource, 0);
    addToRunTimeSelectionTable(massSourceTermModel,cavitationMassSource, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cavitationMassSource::cavitationMassSource
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
    cutoff_(modelDict().lookupOrDefault<scalar>("cutoff",1e-3)),
    n_(modelDict().get<dimensionedScalar>("n")),
    dNuc_(modelDict().get<dimensionedScalar>("dNuc")),
    Cc_(modelDict().getOrDefault<dimensionedScalar>("Cc", dimensionedScalar(dimless, 1))),
    Cv_(modelDict().getOrDefault<dimensionedScalar>("Cv", dimensionedScalar(dimless, 1))),
    p0_(dimensionedScalar("p0", dimPressure, 0))
{}

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::dimensionedScalar Foam::cavitationMassSource::alphaNuc() const
{
    dimensionedScalar Vnuc = n_ * constant::mathematical::pi * pow3(dNuc_) / 6;
    return Vnuc / (1 + Vnuc);
}

Foam::tmp<Foam::volScalarField> Foam::cavitationMassSource::rRb
(
    const volScalarField& limitedAlpha1
) const
{
    tmp<volScalarField> t_rRb = pow
    (
        ((4 * constant::mathematical::pi * n_) / 3)
        * limitedAlpha1 / (1.0 + alphaNuc() - limitedAlpha1),
        1.0/3.0
    );
    
    // 恢复正确量纲：1/m
    t_rRb.ref().dimensions().reset(dimless/dimLength);
    
    return t_rRb;
}

Foam::tmp<Foam::volScalarField> Foam::cavitationMassSource::pCoeff
(
    const volScalarField& p
) const
{
    const volScalarField& rho1 = phase1_.thermo().rho();
    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& pSat = satModel_.pSat();
    const volScalarField& alpha1 = phase1_;

    volScalarField limitedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));
    
    volScalarField rhoMixture = limitedAlpha1 * rho1 + (1 - limitedAlpha1) * rho2;
    
    return
        (3 * rho1 * rho2) * sqrt(2.0/(3.0 * rho1))
        * rRb(limitedAlpha1) / (rhoMixture * sqrt(mag(p - pSat) + 0.01 * pSat));
}

Foam::Pair<Foam::tmp<volScalarField>> Foam::cavitationMassSource::mDotAlphal() const
{
    const volScalarField& p = p_;
    const volScalarField& pSat = satModel_.pSat();
    const volScalarField& alpha1 = phase1_;
    
    volScalarField limitedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));
    //volScalarField limitedAlpha1(clamp(alpha1, zero_one{}));
    volScalarField pCoeffVal = this->pCoeff(p);
    
    return Pair<tmp<volScalarField>>
    (
        Cc_ * limitedAlpha1 * pCoeffVal * max(p - pSat, p0_),
        Cv_ * (1.0 + alphaNuc() - limitedAlpha1) 
            * pCoeffVal * min(p - pSat, p0_)
    );
}

Foam::Pair<Foam::tmp<volScalarField>> Foam::cavitationMassSource::mDotP() const
{
    const volScalarField& p = p_;
    const volScalarField& pSat = satModel_.pSat();
    const volScalarField& alpha1 = phase1_;

    volScalarField pCoeffVal = this->pCoeff(p);
    volScalarField limitedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));
    volScalarField apCoeff = limitedAlpha1 * pCoeffVal;

    return Pair<tmp<volScalarField>>
    (
        Cc_*(1.0 - limitedAlpha1)*pos0(p - pSat)*apCoeff,
        (-Cv_)*(1.0 + alphaNuc() - limitedAlpha1)*neg(p - pSat)*apCoeff
    );
}

// * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * //

Foam::Pair<Foam::tmp<volScalarField>> Foam::cavitationMassSource::vDotAlphal() const
{
    // const volScalarField& p     = p_;
    // const volScalarField& pSat  = satModel_.pSat();
    const volScalarField& rho1  = phase1_.thermo().rho();
    const volScalarField& rho2  = phase2_.thermo().rho();
    const volScalarField& alpha1 = phase1_;

    //volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
    //volScalarField limitedAlpha1(clamp(alpha1, zero_one{}));

    Pair<tmp<volScalarField>> mdot = this->mDotAlphal();

    volScalarField coeff(1.0/rho1 - alpha1*(1.0/rho1 - 1.0/rho2));

    //Info<< "=== CavitationSchnerrSauer::vDotAlphal() Debug ===" << endl;
    //Info<< "mDot dimensions: "  << mdot[0]().dimensions() << endl;
    //Info<< "coeff dimensions: " << coeff.dimensions()     << endl;
    //Info<< "rRb dimensions: "   << rRb(limitedAlpha1)().dimensions() << endl;
    //Info<< "pCoeff dimensions: "<< pCoeff(p, pSat, limitedAlpha1)().dimensions() << endl;

    // Compute raw condensation and evaporation fields
    volScalarField vDotc_field(coeff * mdot[0]());
    volScalarField vDotv_field(coeff * mdot[1]());

    auto makeField = [&](const volScalarField& src)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    src.name(),
                    phase1_.mesh().time().timeName(),
                    phase1_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                src
            )
        );
    };

    return Pair<tmp<volScalarField>>(makeField(vDotc_field), makeField(vDotv_field));
}


Foam::Pair<Foam::tmp<volScalarField>> Foam::cavitationMassSource::vDotP() const
{
    // const volScalarField& p      = p_;
    // const volScalarField& pSat   = satModel_.pSat();
    //const volScalarField& alpha1 = phase1_;
    const volScalarField& rho1   = phase1_.thermo().rho();
    const volScalarField& rho2   = phase2_.thermo().rho();

   // volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    Pair<tmp<volScalarField>> mdot = this->mDotP();

    volScalarField coeff(1.0/rho1 - 1.0/rho2);

    //Info<< "=== CavitationSchnerrSauer::vDotP() Debug ===" << endl;
    //Info<< "mDot dimensions: "  << mdot[0]().dimensions() << endl;
    //Info<< "coeff dimensions: " << coeff.dimensions()     << endl;
    //Info<< "rRb dimensions: "   << rRb(limitedAlpha1)().dimensions() << endl;
    //Info<< "pCoeff dimensions: "<< pCoeff(p, pSat, limitedAlpha1)().dimensions() << endl;

    volScalarField vDotP0_field(coeff * mdot[0]());
    volScalarField vDotP1_field(coeff * mdot[1]());

    auto makeField = [&](const volScalarField& src)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    src.name(),
                    phase1_.mesh().time().timeName(),
                    phase1_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                src
            )
        );
    };

    return Pair<tmp<volScalarField>>(makeField(vDotP0_field), makeField(vDotP1_field));
}

Foam::tmp<Foam::volScalarField> Foam::cavitationMassSource::massSource(volScalarField& rhoSource)
{
    return rhoSource * 0.0;
}

Foam::tmp<Foam::volScalarField> Foam::cavitationMassSource::alphaSource(volScalarField& rhoSource)
{
    return rhoSource / phase1_.thermo().rho() * 0.0;
}

// ************************************************************************* //

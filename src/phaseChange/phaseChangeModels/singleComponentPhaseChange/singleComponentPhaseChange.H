/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Class
    Foam::singleComponentPhaseChange

Description

SourceFiles
    singleComponentPhaseChangeI.H
    singleComponentPhaseChange.C
    singleComponentPhaseChangeIO.C

\*---------------------------------------------------------------------------*/

#ifndef singleComponentPhaseChange_H
#define singleComponentPhaseChange_H

#include "energySourceTermModel.H"
#include "massSourceTermModel.H"
#include "singleComponentSatProp.H"

#include "phaseModel.H"
#include "macroModel.H"
#include "twoPhaseModelThermo.H"
#include "compressibleInterPhaseTransportModel.H"

#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class singleComponentPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class singleComponentPhaseChange
  :   public IOdictionary // will write the pointer in registry
{
    // Private data
    //- Reference to mesh
    const fvMesh& mesh_;

    autoPtr<energySourceTermModel> energyModel_;

    //- Stores the functions
    PtrList<macroModel> macroModels_;

    autoPtr<massSourceTermModel> massModel_;

    autoPtr<singleComponentSatProp> satProp_;

    //- Thermo-package of phase 1
    const phaseModel& phase1_;

    //- Thermo-package of phase 2
    const phaseModel& phase2_;

    // reference to pressure
    const volScalarField& p_;

    const compressibleInterPhaseTransportModel& turbModel_;

    reconstructionSchemes& surf_;

    volScalarField psi0_;

    volScalarField massSource_;

    volScalarField alphaSource_;

    volScalarField vDotAlphaCon_;

    volScalarField vDotAlphaEvp_;

    bool limitHeatFlux_;


public:
    reconstructionSchemes& surf()
    {
        return surf_;
    }

    //- Runtime type information
   // TypeName("singleComponentPhaseChange");

    // Constructors

        // - construct from twoPhaseModelThermo
        singleComponentPhaseChange
        (
            const twoPhaseModelThermo& mixture,
            const volScalarField& p,
            const compressibleInterPhaseTransportModel& turbModel,
            reconstructionSchemes &surf
        );

    //- Destructor
    virtual ~singleComponentPhaseChange();


    // Member Functions


    singleComponentSatProp& satProp()
    {
        return satProp_.ref();
    }

    const phaseModel& phase1()
    {
        return phase1_;
    }

    const phaseModel& phase2()
    {
        return phase2_;
    }

    const compressibleInterPhaseTransportModel& turb()
    {
        return turbModel_;
    }

    void correct();

    // void correct( const surfaceScalarField& kappa1,const surfaceScalarField& kappa2);
    //- Update properties
    void correctSatProperties
    (
        const volScalarField& p,
        const volScalarField& T
    );

    volScalarField& massSource()
    {
        return massSource_;
    }

    volScalarField& alphaSource()
    {
        return alphaSource_;
    }

    //- Return the volumetric condensation and vaporisation rates as a
    //  coefficient to multiply (1 - alphal) for the condensation rate
    //  and a coefficient to multiply  alphal for the vaporisation rate
    volScalarField& vDotAlphaCon()
    {
        return vDotAlphaCon_;
    }
    volScalarField& vDotAlphaEvp()
    {
        return vDotAlphaEvp_;
    }
    //Pair<tmp<volScalarField>> vDotAlphal() const
    //{
    //    return massModel_->vDotAlphal();
    //}

    //- Return the volumetric condensation and vaporisation rates as
    //  coefficients to multiply (p - pSat)
    Pair<tmp<volScalarField>> vDotP() const
    {
        return massModel_->vDotP();
    }

    const volScalarField& psi0()
    {
        return psi0_;
    }

    tmp<fvScalarMatrix> heatTransfer(const volScalarField& T)
    {
        return energyModel_->heatTransfer(T);
    }

    tmp<volScalarField> volTransfer(const volScalarField& T)
    {
        return massModel_->volTransfer(T);
    }

    //tmp<volScalarField> alphaTransfer(const volScalarField& T)
    //{
    //    return massModel_->alphaTransfer(T);
    //}

    tmp<volScalarField> dmdt(const volScalarField& T)
    {
        return massModel_->dmdt(T);
    }

    tmp<fvScalarMatrix> TSource1()
    {
        tmp<fvScalarMatrix> tT1Eqn = energyModel_->TSource1();
        fvScalarMatrix& T1Eqn = tT1Eqn.ref();
        for (auto& mModel: macroModels_)
        {
            mModel.TSource1(T1Eqn);
        }
        return tT1Eqn;
    }

    tmp<fvScalarMatrix> TSource2()
    {
        tmp<fvScalarMatrix> tT2Eqn = energyModel_->TSource2();
        fvScalarMatrix& T2Eqn = tT2Eqn.ref();
        for (auto& mModel: macroModels_)
        {
            mModel.TSource2(T2Eqn);
        }
        return tT2Eqn;
    }

    tmp<volScalarField> energyFlux1()
    {
        return energyModel_->energyFlux1();
    }

    tmp<volScalarField> energyFlux2()
    {
        return energyModel_->energyFlux2();
    }

    tmp<surfaceScalarField> getalphaf()
    {
        return energyModel_->getalphaf();
    }



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

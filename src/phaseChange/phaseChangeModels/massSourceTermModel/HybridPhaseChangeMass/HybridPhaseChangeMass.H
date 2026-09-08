/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR),
                                     Univeristy College Dublin(UCD)
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

Class
    Foam::HybridPhaseChangeMass

Description
    hybrid phase change mass source model 

SourceFiles
    HybridPhaseChangeMass.C

\*---------------------------------------------------------------------------*/

#ifndef HybridPhaseChangeMass_H
#define HybridPhaseChangeMass_H

#include "typeInfo.H"
#include "massSourceTermModel.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "fvScalarMatrix.H"
#include "fvm.H"

#include "reconstructionSchemes.H"
#include "singleComponentSatProp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class HybridPhaseChangeMass Declaration
\*---------------------------------------------------------------------------*/

class HybridPhaseChangeMass
:
    public massSourceTermModel
{
    // Private data
        
        //- Sub-model for interface region (e.g., Schrage, implicitGrad)
        autoPtr<massSourceTermModel> interfaceModel_;
        
        //- Sub-model for liquid region (cavitation model)
        autoPtr<massSourceTermModel> cavitationModel_;


    // Private Member Functions
        //- Disallow copy construct
        HybridPhaseChangeMass(const HybridPhaseChangeMass&);
        
        //- Disallow default bitwise assignment
        void operator=(const HybridPhaseChangeMass&);

        //- Zero out interface-region contributions in cavitation source
        //  Prevents double-counting with interfaceModel in interface cells
        void zeroInterfaceCells(volScalarField& field) const;

public:
    //- Runtime type information
    TypeName("HybridPhaseChangeMass");


    // Constructors
        //- Construct from components
        HybridPhaseChangeMass
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );


    //- Destructor
    virtual ~HybridPhaseChangeMass()
    {}


    // Member Functions   
    virtual tmp<volScalarField> massSource(volScalarField& rhoSource);

    virtual tmp<volScalarField> alphaSource(volScalarField& rhoSource);

    //- Return the mass condensation and vaporisation rates as a
    //  coefficient to multiply (1 - alphal) for the condensation rate
    //  and a coefficient to multiply  alphal for the vaporisation rate
    virtual Pair<tmp<volScalarField>> vDotAlphal() const;

    //- Return the mass condensation and vaporisation rates as coefficients
    //  to multiply (p - pSat)
    virtual Pair<tmp<volScalarField>> vDotP() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

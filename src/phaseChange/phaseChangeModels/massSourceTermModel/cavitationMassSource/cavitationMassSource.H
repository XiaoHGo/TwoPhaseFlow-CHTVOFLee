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

Class
    Foam::cavitationMassSource

Description
    cavitationMassSource mass source terms are zero

SourceFiles
    cavitationMassSource.C

\*---------------------------------------------------------------------------*/

#ifndef cavitationMassSource_H
#define cavitationMassSource_H

#include "typeInfo.H"
#include "massSourceTermModel.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "fvScalarMatrix.H"
#include "fvm.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class cavitationMassSource Declaration
\*---------------------------------------------------------------------------*/

class cavitationMassSource
    : public massSourceTermModel
{
private:
    // Private data
        scalar cutoff_;
        // SchnerrSauer 模型参数
        //- Bubble number density (nuclei concentration)
        dimensionedScalar n_;
        
        //- Nucleation site diameter
        dimensionedScalar dNuc_;
        
        //- Condensation rate coefficient
        dimensionedScalar Cc_;
        
        //- Vaporisation rate coefficient
        dimensionedScalar Cv_;
        
        //- Reference pressure (for model compatibility)
        dimensionedScalar p0_;
        
        //- Specific gas constant (for incompressible case)
        dimensionedScalar R_;

    // Private Member Functions

        //- Disallow copy construct
        cavitationMassSource(const cavitationMassSource&);

        //- Disallow default bitwise assignment
        void operator=(const cavitationMassSource&);
        
        // 私有辅助函数
        dimensionedScalar alphaNuc() const;
        tmp<volScalarField> rRb(const volScalarField& limitedAlpha1) const;
        tmp<volScalarField> pCoeff(const volScalarField& p) const;
    
        // 计算质量相变率 (kg/m²·s)
        Pair<tmp<volScalarField>> mDotAlphal() const;
        Pair<tmp<volScalarField>> mDotP() const;

public:

    //- Runtime type information
    TypeName("cavitationMassSource");

    // Constructors

        //- Construct from components
        cavitationMassSource
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    //- Destructor
    virtual ~cavitationMassSource()
    {}

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

/*---------------------------------------------------------------------------*\
            Copyright (c) 2025-2028, University College Dublin (UCD)

            Extended Lee boiling mass-source model
\*---------------------------------------------------------------------------*/

#ifndef boilingLeeMassSource_H
#define boilingLeeMassSource_H

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
                  Class boilingLeeMassSource Declaration
\*---------------------------------------------------------------------------*/

class boilingLeeMassSource
:
    public massSourceTermModel
{
private:

    // Private data
        //- Lee relaxation coefficient [1/s]
        dimensionedScalar C_;
        //- Activation temperature [K]
        dimensionedScalar Tactivate_;
        //- Minimum donor phase volume fraction
        scalar alphaMin_;

    // Private Member Functions
        //- Explicit Lee volumetric mass-transfer rate [kg/m3/s]
        //  phase1 = liquid
        //  mDot = C * alpha1 * rho1 * (T - Tactivate)/Tactivate
        //  for T > Tactivate when C > 0.
        tmp<volScalarField> Kexp(const volScalarField& T) const;
        //- Disallow copy construct
        boilingLeeMassSource(const boilingLeeMassSource&);
        //- Disallow default bitwise assignment
        void operator=(const boilingLeeMassSource&);

public:

    //- Runtime type information
    TypeName("boilingLeeMassSource");

    // Constructors
        boilingLeeMassSource
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    //- Destructor
    virtual ~boilingLeeMassSource()
    {}

    // Member Functions
        // Lee-specific interfaces
        virtual tmp<volScalarField> volTransfer(const volScalarField& T) const;
        //virtual tmp<volScalarField> alphaTransfer(const volScalarField& T) const;
        virtual tmp<volScalarField> dmdt(const volScalarField& T) const;

        virtual tmp<volScalarField> massSource(volScalarField& rhoSource);
        virtual tmp<volScalarField> alphaSource(volScalarField& rhoSource);

        //virtual Pair<tmp<volScalarField>> vDotAlphal() const;
        //virtual Pair<tmp<volScalarField>> vDotP() const;

};


} // End namespace Foam

#endif

// ************************************************************************* //

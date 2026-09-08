
#ifndef BoilingLee_H
#define BoilingLee_H

#include "typeInfo.H"
#include "energySourceTermModel.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "fvScalarMatrix.H"
#include "fvm.H"

#include "reconstructionSchemes.H"
#include "singleComponentSatProp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class BoilingLee Declaration
\*---------------------------------------------------------------------------*/

class BoilingLee
:
    public energySourceTermModel
{
    // Private data
        //- Lee relaxation coefficient [1/s]
        dimensionedScalar C_;

        //- Activation temperature [K]
        dimensionedScalar Tactivate_;

        //- Minimum liquid fraction required for phase change
        scalar alphaMin_;

    // Private Member Functions
        //- Full explicit Lee mass-transfer rate [kg/m3/s]
        tmp<volScalarField> Kexp(const volScalarField& T) const;

        //- Temperature-linearised coefficient
        //  mDot = KSp*T + KSu
        //  Units [kg/m3/s/K]
        tmp<volScalarField> KSp(const volScalarField& T) const;

        //- Explicit part of temperature-linearised mass transfer
        //  Units [kg/m3/s]
        tmp<volScalarField> KSu(const volScalarField& T) const;

        //- Latent heat used by Lee energy coupling [J/kg]
        tmp<volScalarField> latentHeat(const volScalarField& T) const;

        //- Disallow copy construct
        BoilingLee(const BoilingLee&);
        
        //- Disallow default bitwise assignment
        void operator=(const BoilingLee&);
        

public:
    //- Runtime type information
    TypeName("BoilingLee");


    // Constructors
        //- Construct from components
        BoilingLee
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const compressibleInterPhaseTransportModel& turbModel,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BoilingLee()
    {}


    // Member Functions
        //- Semi-implicit latent-heat source for global T equation
        virtual tmp<fvScalarMatrix> heatTransfer(const volScalarField& T);

        //- Return energy source term for phase 1 temperature equation
        virtual tmp<fvScalarMatrix> TSource1();
        
        //- Return energy source term for phase 2 temperature equation
        virtual tmp<fvScalarMatrix> TSource2();
        
        //- Return total energy source (q_pc in W/m³)
        virtual tmp<volScalarField> energySource();
        
        //- Return energy flux for phase 1 (not used for cavitation)
        virtual tmp<volScalarField> energyFlux1();
        
        //- Return energy flux for phase 2 (not used for cavitation)
        virtual tmp<volScalarField> energyFlux2();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

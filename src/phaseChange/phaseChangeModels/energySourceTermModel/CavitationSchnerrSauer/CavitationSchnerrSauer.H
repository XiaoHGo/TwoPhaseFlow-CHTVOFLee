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
    Foam::Schrage

Description
    interface heat resistance model based on the Schrage equation

SourceFiles
    Schrage.C

\*---------------------------------------------------------------------------*/
#ifndef CavitationSchnerrSauer_H
#define CavitationSchnerrSauer_H

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
                Class CavitationSchnerrSauer Declaration
\*---------------------------------------------------------------------------*/

class CavitationSchnerrSauer
:
    public energySourceTermModel
{
    // Private data
    // Private Member Functions
        //- Disallow copy construct
        CavitationSchnerrSauer(const CavitationSchnerrSauer&);
        
        //- Disallow default bitwise assignment
        void operator=(const CavitationSchnerrSauer&);
        

public:
    //- Runtime type information
    TypeName("CavitationSchnerrSauer");


    // Constructors
        //- Construct from components
        CavitationSchnerrSauer
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
    virtual ~CavitationSchnerrSauer()
    {}


    // Member Functions
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

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

#include "massSourceTermModel.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(massSourceTermModel, 0);
    defineRunTimeSelectionTable(massSourceTermModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massSourceTermModel::massSourceTermModel
(
    const word& type,
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    reconstructionSchemes& surf,
    const dictionary& dict
)
:
    dictionary(dict),
    massSourceTermModelCoeffs_(dict.optionalSubDict(type + "Coeffs")),
    phase1_(phase1),
    phase2_(phase2),
    p_(p),
    satModel_(satModel),
    surf_(surf)
{

}


const Foam::dictionary&
Foam::massSourceTermModel::modelDict() const
{
    return massSourceTermModelCoeffs_;
}


Foam::dictionary&
Foam::massSourceTermModel::modelDict()
{
    return massSourceTermModelCoeffs_;
}

Foam::tmp<Foam::volScalarField>Foam::massSourceTermModel::dmdt(const volScalarField& T) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "zeroDmdt",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase1_.mesh(),
            dimensionedScalar
            (
                dimDensity/dimTime,
                Zero
            )
        )
    );
}

Foam::tmp<Foam::volScalarField>Foam::massSourceTermModel::volTransfer(const volScalarField& T) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "zeroVolTransfer",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase1_.mesh(),
            dimensionedScalar
            (
                dimless/dimTime,
                Zero
            )
        )
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>Foam::massSourceTermModel::vDotAlphal() const
{
    auto makeZero = [&]()
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "zeroVdotAlphal",
                    phase1_.mesh().time().timeName(),
                    phase1_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phase1_.mesh(),
                dimensionedScalar(dimless/dimTime, 0)
            )
        );
    };

    return Pair<tmp<volScalarField>>(makeZero(), makeZero());
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>Foam::massSourceTermModel::vDotP() const
{
    auto makeZero = [&]()
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "zeroVdotP",
                    phase1_.mesh().time().timeName(),
                    phase1_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phase1_.mesh(),
                dimensionedScalar
                (
                    dimTime/(dimDensity*dimLength*dimLength),  // [-1 1 1 0 0 0 0]
                    0
                )
            )
        );
    };

    return Pair<tmp<volScalarField>>(makeZero(), makeZero());
}


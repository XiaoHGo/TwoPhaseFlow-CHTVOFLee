volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);

const volScalarField& alphaSource = phaseChange.alphaSource();

forAll(dgdt, celli)
{
    if (dgdt[celli] > 0.0)
    {
        Sp[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
        Su[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
    }
    else if (dgdt[celli] < 0.0)
    {
        Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
    }
    if (alphaSource[celli] > 0.0 && alpha1[celli] > 0.0)
    {
        //Sp[celli] -= alphaSource[celli];
        Su[celli] -= alphaSource[celli];
    }
    else if (alphaSource[celli] < 0.0 && alpha1[celli] < 1.0)
    {
        Su[celli] -= alphaSource[celli];
    }

}

// ========== add: the cavitation source term ==========
//Pair<tmp<volScalarField>> vDotAlphal = phaseChange.vDotAlphal();
volScalarField& vDotc = phaseChange.vDotAlphaCon();
volScalarField& vDotv = phaseChange.vDotAlphaEvp();

// transfer internal class
volScalarField::Internal Sp_cav = (vDotv - vDotc)();
volScalarField::Internal Su_cav = vDotc();
// ========== add: the cavitation source term =======
forAll(Sp_cav, celli)
{
    Sp[celli] += Sp_cav[celli];
    Su[celli] += Su_cav[celli];
}

volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phi + mesh.phi())
  : fvc::div(phi)
);

// ========== add: the Lee source term ==========
tmp<volScalarField> tdmdtLee = phaseChange.dmdt(T);

const volScalarField& dmdtNet = tdmdtLee();

// ========== add: the Lee source term =======
#include "alphaSuSpArrange.H"

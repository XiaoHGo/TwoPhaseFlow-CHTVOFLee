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

// ========== 新增：获取空化源项 ==========
//Pair<tmp<volScalarField>> vDotAlphal = phaseChange.vDotAlphal();
volScalarField& vDotc = phaseChange.vDotAlphaCon();
volScalarField& vDotv = phaseChange.vDotAlphaEvp();

// 转换 Internal 类型
volScalarField::Internal Sp_cav = (vDotv - vDotc)();
volScalarField::Internal Su_cav = vDotc();
// ========== 叠加空化源项 ==========
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

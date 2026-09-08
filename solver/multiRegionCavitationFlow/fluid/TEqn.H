// ==================================================================
// TEqns.H
// ==================================================================

// ------------------------------------------------------------------
// 1. Thermal properties
// ------------------------------------------------------------------
// rhoCp  = alpha1 * rhoCpL + alpha2 * rhoCpV;
rhoCp  = rhoCpL + rhoCpV;

surfaceScalarField kappaEff
(
    "kappaEff",
    fvc::interpolate(alpha1)
   *fvc::interpolate(mixture.thermo1().kappaEff(turbulence.alphat()))
  + fvc::interpolate(alpha2)
   *fvc::interpolate(mixture.thermo2().kappaEff(turbulence.alphat()))
);

// rhoCpPhi = fvc::interpolate(alpha1) * rhoCpLPhi + fvc::interpolate(alpha2) * rhoCpVPhi;
rhoCpPhi =  rhoCpLPhi + rhoCpVPhi;

const volScalarField& TSat = phaseChange.satProp().TSat();
const volScalarField& pSat = phaseChange.satProp().pSat();
//const volScalarField& L    = phaseChange.satProp().L();

Info<< "=== TEqns Debug ===" << endl;
Info<< "rhoCpL min/max: " << gMin(rhoCpL) << " / " << gMax(rhoCpL) << endl;
Info<< "rhoCpV min/max: " << gMin(rhoCpV) << " / " << gMax(rhoCpV) << endl;
Info<< "TSat min/max:   " << gMin(TSat)   << " / " << gMax(TSat)   << endl;
Info<< "pSat min/max:   " << gMin(pSat)   << " / " << gMax(pSat)   << endl;
Info<< "massSource min/max: "
    << gMin(phaseChange.massSource()) << " / "
    << gMax(phaseChange.massSource()) << endl;

// ------------------------------------------------------------------
// 2. Global mixture temperature equation
// ------------------------------------------------------------------
fvScalarMatrix TEqn
(
    fvm::ddt(rhoCp, T)
  + fvm::div(rhoCpPhi, T)
  - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
  - fvm::laplacian(kappaEff, T)
 ==
    fvOptions(rhoCp, T)
);

TEqn.relax();
TEqn.solve();

T = max(min(T, TMax), TMin);
T.correctBoundaryConditions();

Info<< "min(T)  " << gMin(T)  << "  max(T)  " << gMax(T)  << endl;

TL=T;
TV=T;

TL.correctBoundaryConditions();
TV.correctBoundaryConditions();

//TL = max(min(TL,TMax),TMin);
//TV = max(min(TV,TMax),TMin);

Info << "min(TL) " << gMin(TL) << "  max(TL) " << gMax(TL)  << endl;
Info << "min(TV) " << gMin(TV) << "  max(TV) " << gMax(TV) << endl;

T.correctBoundaryConditions();

mixture.correctThermo(TL,TV);

mixture.correct();

//#include "TvlEqns.H"
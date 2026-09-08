// ------------------------------------------------------------------
// 3. Interface reconstruction: TV from gas-side heat flux,
//    TL from energy conservation
//
// Gas phase (Fo >> 1, quasi-steady linear):
//   qV = energyFlux2() [W/m²], already computed by erfStefan
//   TV_avg = TSat + qV * d_V / (2*kV)
//
// Liquid phase (energy conservation):
//   TL_avg = (T_mix - TV_avg * alpha2) / alpha1
//
// Pure-phase cells: TL = TV = T (no reconstruction needed)
// ------------------------------------------------------------------
scalar cutOff = 1e-3;

// Get gas-side heat flux from erfStefan (quasi-steady linear interp)
volScalarField qV_field = phaseChange.energyFlux2();
//const volScalarField& TSat = phaseChange.satProp().TSat();
//const volScalarField& L    = phaseChange.satProp().L();
reconstructionSchemes& surf = phaseChange.surf();

forAll(alpha1, cellI)
{
    const scalar a = alpha1[cellI];

    if (!surf.interfaceCell()[cellI] || a < cutOff || a > (1-cutOff))
    {
        // Pure-phase or non-interface: follow global T
        TL[cellI] = T[cellI];
        TV[cellI] = T[cellI];
        continue;
    }

    // Interface cell geometry
    const scalar V_cell = mesh.V()[cellI];
    const scalar A_s    = mag(surf.normal()[cellI]);
    if (A_s < SMALL)
    {
        TL[cellI] = T[cellI];
        TV[cellI] = T[cellI];
        continue;
    }

    const scalar d_V  = (1.0-a) * V_cell / A_s;
    const scalar Ts   = TSat[cellI];
    const scalar Tmix = T[cellI];
    const scalar kV   = mixture.thermo2().kappa()()[cellI];
    const scalar qV   = qV_field[cellI];  // [W/m²]

    // 检查温差是否足够大
    // 如果T接近TSat，界面处已经接近热平衡，不需要重建
    scalar dT = mag(Tmix - Ts);
    if (dT < 0.5)  // 温差小于0.5K，认为无驱动力
    {
        TL[cellI] = Tmix;
        TV[cellI] = Tmix;
        continue;
    }// 有温差，才做界面重建

    // Gas side: quasi-steady linear profile
    // TV_avg = TSat + qV * (d_V/2) / kV
    // d_V/2 is the distance from interface to gas sub-volume centroid
    const scalar TV_avg = Ts + qV * d_V / (2.0*kV + SMALL);

    // Liquid side: energy conservation
    // TL_avg * alpha1 + TV_avg * alpha2 = T_mix
    const scalar TL_avg = (Tmix - TV_avg*(1.0-a)) / (a + SMALL);

    if (TL_avg < TMin.value() || TL_avg > TMax.value() || !std::isfinite(TL_avg))
    {
        TL[cellI] = Tmix;
        TV[cellI] = Tmix;
        continue;
    }

    TL[cellI] = max(min(TL_avg, TMax.value()), TMin.value());
    TV[cellI] = max(min(TV_avg, TMax.value()), TMin.value());
}

TL.correctBoundaryConditions();
TV.correctBoundaryConditions();

// ------------------------------------------------------------------
// 4. Update global T from energy conservation
// T = (rhoCpL*TL*alpha1 + rhoCpV*TV*alpha2) / (rhoCpL*alpha1 + rhoCpV*alpha2)
// ------------------------------------------------------------------
T = ( rhoCpL*TL*alpha1 + rhoCpV*TV*alpha2 )
  / max
    (
        rhoCpL*alpha1 + rhoCpV*alpha2,
        dimensionedScalar("rhoCpMin", rhoCp.dimensions(), SMALL)
    );

T = max(min(T, TMax), TMin);
T.correctBoundaryConditions();

Info<< "min(TL) " << gMin(TL) << "  max(TL) " << gMax(TL) << endl;
Info<< "min(TV) " << gMin(TV) << "  max(TV) " << gMax(TV) << endl;
Info<< "min(T)  " << gMin(T)  << "  max(T)  " << gMax(T)  << endl;

// ------------------------------------------------------------------
// 5. Update thermophysical properties
// ------------------------------------------------------------------
mixture.correctThermo(TL, TV);
mixture.correct();
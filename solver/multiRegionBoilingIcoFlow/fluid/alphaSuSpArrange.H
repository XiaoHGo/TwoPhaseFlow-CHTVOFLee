// ==================================================================
// alphaSuSpArrange.H
//
// Arrange phase-change mass transfer into Su/Sp for phase1
//
// Sign convention:
//
//      dmdtNet > 0 : mass transfer INTO phase1
//      dmdtNet < 0 : mass transfer OUT OF phase1
//
// For phase1 = liquid:
//
//      condensation : dmdtNet > 0
//      evaporation  : dmdtNet < 0
//
// Based on icoReactingMultiphaseInterFoam /
// MassTransferPhaseSystem::alphaTransfer()
// ==================================================================

const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();

volScalarField coeff1("alphaTransferCoeff1", 1.0/rho1);
volScalarField coeff2("alphaTransferCoeff2",1.0/rho2);
volScalarField coeff12("alphaTransferCoeff12",coeff1 - coeff2);

forAll(dmdtNet, celli)
{
    const scalar mdot = dmdtNet[celli];

    const scalar c1  = coeff1[celli];
    const scalar c12 = coeff12[celli];

    const scalar alpha1Limited = min(max(alpha1[celli], scalar(0)), scalar(1));

    // --------------------------------------------------------------
    // Direct phase-mass source
    //
    //      dmdtNet/rho1
    // --------------------------------------------------------------

    Su[celli] += c1*mdot;

    // --------------------------------------------------------------
    // Volume-change correction
    //
    //     -alpha1*dmdtNet*(1/rho1 - 1/rho2)
    //
    // Distributed into Su/Sp according to sign for boundedness.
    // --------------------------------------------------------------

    if (mdot > 0.0)
    {
        if (c12 > 0.0)
        {
            // Implicit treatment
            Sp[celli] -= mdot*c12;
        }
        else if (c12 < 0.0)
        {
            // Explicit treatment
            Su[celli] -= mdot*c12*alpha1Limited;
        }
    }
    else if (mdot < 0.0)
    {
        if (c12 > 0.0)
        {
            // Explicit treatment
            Su[celli] -= mdot*c12*alpha1Limited;
        }
        else if (c12 < 0.0)
        {
            // Implicit treatment
            Sp[celli] -= mdot*c12;
        }
    }
}

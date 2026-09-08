Info<< nl<< "========== Global continuity debug =========="<< endl;
// --------------------------------------------------
// 1. Total mixture mass
// --------------------------------------------------
//const volScalarField& rho1 = mixture.thermo1().rho();
//const volScalarField& rho2 = mixture.thermo2().rho();
const scalar M1 = gSum(alpha1.primitiveField()*rho1.primitiveField()*mesh.V());
const scalar M2 = gSum(alpha2.primitiveField()*rho2.primitiveField()*mesh.V());
const scalar MtotNow = M1 + M2;
// --------------------------------------------------
// 2. Lee phase change source term
// --------------------------------------------------
tmp<volScalarField> tdmdtMass = phaseChange.dmdt(T);
const volScalarField& dmdtMass = tdmdtMass();
    //scalar M1 = gSum(alpha1.primitiveField()*mixture.thermo1().rho().primitiveField()*mesh.V());
    //scalar M2 = gSum(alpha2.primitiveField()*mixture.thermo2().rho().primitiveField()*mesh.V());
scalar mdotInt =gSum(dmdtMass.primitiveField()*mesh.V());
    //Info<< "M1 = " << M1
    //    << ", M2 = " << M2
    //    << ", Mtot = " << M1 + M2
    //    << ", integral(dmdt dV) = " << mdotInt
    //    << endl;
// --------------------------------------------------
// 3. External-boundary mass flux
// --------------------------------------------------
scalar boundaryMassFlux = 0.0;
forAll(rhoPhi.boundaryField(), patchi)
{
    const fvPatch& patch = mesh.boundary()[patchi];
    // Skip processor/cyclic and other coupled internal boundaries
    if (!patch.coupled())
    {
        boundaryMassFlux +=
            sum(rhoPhi.boundaryField()[patchi]);
    }
}
// Exactly one MPI collective call
reduce(boundaryMassFlux, sumOp<scalar>());
// --------------------------------------------------
// 4. Global continuity residual
// --------------------------------------------------
const scalar dt = mesh.time().deltaTValue();
const scalar continuityResidual = mdotInt + boundaryMassFlux;
const scalar massScale = MtotNow/dt;
// --------------------------------------------------
// 5. Print
// --------------------------------------------------
Info<< "M1= "<< M1 << " kg" << nl
    << "M2= "<< M2 << " kg" << nl
    << "Mtot= "<< MtotNow << " kg" << nl
    << "integral(dmdt dV)= "<< mdotInt << " kg/s/" << nl
    << "Boundary int(rhoPhi)= "<< boundaryMassFlux << " kg/s" << nl
    << "Continuity residual= "<< continuityResidual << " kg/s" << nl
    << "Normalized residual= "<< continuityResidual/(massScale + SMALL)
    << endl;
// Update for next timestep
// previousFluidMass = MtotNow;

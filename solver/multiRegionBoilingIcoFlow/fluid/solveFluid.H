if (finalIter)
{
    //mesh.data::add("finalIteration", true);
    //for OpenFOAM v2206 or earlier
    mesh.data().setFinalIteration(true); // for OpenFOAM v2406 or later
}


//    #include "readDyMControls.H"
    #include "readTimeControls.H"

    bool correctPhi = pimple.lookupOrDefault
    (
        "correctPhi",
        correctPhi
    );

    // Store divU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
   volScalarField absdivU("divU0", fvc::div(fvc::absolute(phi, U)));

    Info<< "Time = " << runTime.timeName() << nl << endl;

// --- Pressure-velocity PIMPLE corrector loop
    if (oCorr == 0)
    {
        scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();
        advector.surf().reconstruct();

        mesh.update();

        if (mesh.changing())
        {
            MRF.update();

            advector.surf().mapAlphaField();
            alpha2 = 1.0 - alpha1;
            alpha2.correctBoundaryConditions();
            rho == alpha1*rho1 + alpha2*rho2;
            rho.correctBoundaryConditions();
            rho.oldTime() = rho;
            alpha2.oldTime() = alpha2;

            Info<< "Execution time for mesh.update() = "
                << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                << " s" << endl;

            // gets recompute by surfaces forces
            // gh = (g & mesh.C()) - ghRef;
            // ghf = (g & mesh.Cf()) - ghRef;
        }

        if ((mesh.changing() && correctPhi))
        {
            // Calculate absolute flux from the mapped surface velocity
            phi = mesh.Sf() & UfFluid[i]();

            #include "correctPhi.H"

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

            mixture.correct();
        }

        // if (mesh.changing() && checkMeshCourantNo)
        // {
        //     #include "meshCourantNo.H"
        // }
    }

    if(overwrite)
    {
        continue;
    }

    #include "alphaControls.H"
    #include "compressibleAlphaEqnSubCycle.H"

    // ========== Mass Conservation debug ==========
    //Info<< "\n ========== Mass Conservation debug ========== \n"<< endl;
    //tmp<volScalarField> tdmdtMass = phaseChange.dmdt(T);
    //const volScalarField& dmdtMass = tdmdtMass();
    //scalar M1 = gSum(alpha1.primitiveField()*mixture.thermo1().rho().primitiveField()*mesh.V());
    //scalar M2 = gSum(alpha2.primitiveField()*mixture.thermo2().rho().primitiveField()*mesh.V());
    //scalar mdotInt =gSum(dmdtMass.primitiveField()*mesh.V());
    //Info<< "M1 = " << M1
    //    << ", M2 = " << M2
    //    << ", Mtot = " << M1 + M2
    //    << ", integral(dmdt dV) = " << mdotInt
    //    << endl;
    //============================================
   
    mixture.correct();

    turbulence.correctPhasePhi();
    surfForces.correct();

    phaseChange.correctSatProperties(p, T);

    #include "YEqns.H"

    #include "UEqn.H"
    #include "TEqn.H"

    // ========== Heat Conservation debug ==========
    //Info<< "\n========== Heat Conservation debug ==========\n" << endl;
    //tmp<volScalarField> tdmdtHeat = phaseChange.dmdt(T);
    //const volScalarField& dmdtHeat = tdmdtHeat();
    //const volScalarField& L = phaseChange.satProp().L();
    //tmp<volScalarField> tqLatent = -dmdtHeat*L;
    //const volScalarField& qLatent = tqLatent();
    //scalar Qlatent = gSum(qLatent.primitiveField()*mesh.V());
    //Info<< "Lee latent power = "
    //    << Qlatent << " W"
    //    << endl;
    //============================================

    phaseChange.correct();

    //const auto vDot = phaseChange.vDotAlphal();
    //Info<< "vDotAlphal ConMin/ConMax | EvpMin/EvpMax: " << endl;
    //Info<< gMin(phaseChange.vDotAlphaCon()) << " / " << gMax(phaseChange.vDotAlphaCon()) << " | "
    //    << gMin(phaseChange.vDotAlphaEvp()) << " / " << gMax(phaseChange.vDotAlphaEvp()) << endl;

    bool closedVolume = p_rgh.needReference();

    #include "correctMass.H"

    // --- Pressure corrector loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"

        // ==========Volume Change debug ==========
        //tmp<volScalarField> tVol =phaseChange.volTransfer(T);
        //const volScalarField& volTransfer =tVol();
        //scalar volInt = gSum(volTransfer.primitiveField()*mesh.V());
        //Info<< "volTransfer min/max = "
        //    << gMin(volTransfer) << " "
        //    << gMax(volTransfer)
        //    << ", int(volTransfer)="
        //    << volInt
        //    << endl;
        //============================================ 
        
    }

    // ==========Global continuity debug ==========
    #include "ContinuityDebug.H"
    //============================================ 

    tCpl = mixture.thermo1().Cp();
    tCpv = mixture.thermo2().Cp();

    rhoCpL = alpha1*rho1*tCpl + alpha2*rho1*tCpl*1e-7;
    rhoCpV = alpha2*rho2*tCpv + alpha1*rho2*tCpv*1e-7;

//      if (pimple.turbCorr())
    {
        turbulence.correct();
    }

if (finalIter)
{
    //mesh.data::remove("finalIteration");  
    //for OpenFOAM v2206 or earlier
    mesh.data().setFinalIteration(false);   // for OpenFOAM v2406 or later
}

if (finalIter)
{
    //mesh.data::add("finalIteration", true);
    //for OpenFOAM v2206 or earlier
    mesh.data().setFinalIteration(true); // for OpenFOAM v2406 or later
}

if(overwrite)
{
    continue;
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(betav*rho, h)
          - (
               thermo.isotropic()
             ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
             : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
          ==
            fvOptions(rho, h)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve(mesh.solver(h.select(finalIter)));

        fvOptions.correct(h);
    }

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}

if (finalIter)
{
    //mesh.data::remove("finalIteration");  
    //for OpenFOAM v2206 or earlier
    mesh.data().setFinalIteration(false);   // for OpenFOAM v2406 or later
}

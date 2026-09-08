{
    PtrList<volScalarField>& Y1 = mixture.phase1().Y();
    PtrList<volScalarField>& Y2 = mixture.phase2().Y();

    volScalarField muEff(turbulence.muEff());

    // ---------------------------------------------------------------------
    // Phase 1 species
    // ---------------------------------------------------------------------
    if (!Y1.empty())
    {
        label inertIndex1 = mixture.phase1().inertIndex();
        volScalarField Yt(0.0*Y1[0]);

        forAll(Y1, i)
        {
            tmp<fvScalarMatrix> Y1iEqn
            (
                mixture.phase1().YiEqn(Y1[i], muEff)
            );

            if (Y1iEqn.valid())
            {
                Y1iEqn->relax();
                Y1iEqn->solve(mesh.solver("Yi"));

                Y1[i].max(0.0);
                Y1[i].min(1.0);

                Yt += Y1[i];
            }
        }

        Y1[inertIndex1] = scalar(1) - Yt;
        Y1[inertIndex1].max(0.0);
        Y1[inertIndex1].min(1.0);
    }


    // ---------------------------------------------------------------------
    // Phase 2 species: gas = air + vapour
    // ---------------------------------------------------------------------
    if (!Y2.empty())
    {
        label inertIndex2 = mixture.phase2().inertIndex();
        volScalarField Yt(0.0*Y2[0]);

        // Lee convention:
        //
        // dmdt > 0 : mass transfer into phase1
        // dmdt < 0 : evaporation from phase1 into phase2
        //
        // Therefore gas mass source during evaporation is:
        //
        //     Sg = -dmdt > 0
        //
        tmp<volScalarField> tdmdt = phaseChange.dmdt(T);
        const volScalarField& dmdt = tdmdt();

        tmp<volScalarField> tSg = -dmdt;
        const volScalarField& Sg = tSg();


        forAll(Y2, i)
        {
            tmp<fvScalarMatrix> Y2iEqn
            (
                mixture.phase2().YiEqn(Y2[i], muEff)
            );

            if (Y2iEqn.valid())
            {
                // Water transferred from liquid phase enters the gas
                // specifically as the "vapour" component.
                //
                // Conservative mass-fraction source:
                //
                //     S_Yv = Sg*(1 - Yv)
                //          = Sg - Sg*Yv
                //
                // Treat -Sg*Yv implicitly for stability.
                //
                if (Y2[i].member() == "vapour")
                {
                    Y2iEqn.ref() +=
                        Sg
                      - fvm::Sp(Sg, Y2[i]);
                }

                Y2iEqn->relax();
                Y2iEqn->solve(mesh.solver("Yi"));

                Y2[i].max(0.0);
                Y2[i].min(1.0);

                Yt += Y2[i];
            }
        }

        // air is the inert/dependent component:
        //
        //     Yair = 1 - sum(Yactive)
        //
        Y2[inertIndex2] = scalar(1) - Yt;
        Y2[inertIndex2].max(0.0);
        Y2[inertIndex2].min(1.0);
    }
}

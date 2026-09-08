
    scalar cutOff = 1e-3;
    //scalar omega_= 0.5;

    surfaceScalarField kappaLEff
    (
        "kappaLEff",
        fvc::interpolate(mixture.thermo1().kappaEff(turbulence.alphat()))
    );
    surfaceScalarField kappaVEff
    (
        "kappaVEff",
        fvc::interpolate(mixture.thermo2().kappaEff(turbulence.alphat()))
    );
    
    //----------------------------//
    // 在 TLEqn 和 TVEqn 的定义之前
    //Pair<tmp<volScalarField>> vDotAlphal = phaseChange.vDotAlphal();
    volScalarField& vDotc = phaseChange.vDotAlphaCon();
    volScalarField& vDotv = phaseChange.vDotAlphaEvp();
    const volScalarField L = phaseChange.satProp().L();

    // 构建液相能量方程的潜热源项 (W/m³)
    // 蒸发时 (vDotv > 0)，液体失去能量 -> 负源项
    volScalarField latentSourceL 
    = (alpha1 * vDotv + alpha2 * vDotc) 
    * rho1 * rho2 / (alpha1 * rho1 + alpha2 * rho2 ) 
    * L;

    // 构建气相能量方程的潜热源项 (W/m³)
    // 蒸发时 (vDotv > 0)，气体获得能量 -> 正源项
    volScalarField latentSourceV 
    = - (alpha1 * vDotv + alpha2 * vDotc) 
    * rho1 * rho2 / (alpha1 * rho1 + alpha2 * rho2 ) 
    * L;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const volScalarField& TSat = phaseChange.satProp().TSat();

    forAll(kappaLEff,faceI)
    {
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];
        const scalar alphaOwn = alpha1[own];
        const scalar alphaNei = alpha1[nei];

        if(alphaOwn < cutOff || alphaNei < cutOff) // liquid phase
        {
            kappaLEff[faceI] = 0.0;
        }
        if(alphaOwn > (1-cutOff) || alphaNei > (1-cutOff)) // gas phase
        {
            kappaVEff[faceI] = 0.0;
        }
    }

    // energy equation
    /*
    fvScalarMatrix TEqn
	(
		fvm::ddt(T)
		+ fvm::div(phi, T)
		- fvm::laplacian(DT, T)
	);
	TEqn.solve();
    */

    fvScalarMatrix TLEqn
    (
        fvm::ddt(rhoCpL,TL)
        + fvm::div(rhoCpLPhi,TL)
	    - fvm::Sp(fvc::ddt(rhoCpL) + fvc::div(rhoCpLPhi), TL)
        - fvm::laplacian(kappaLEff,TL)
        + phaseChange.TSource1()
        - latentSourceL
            ==
        fvOptions(rhoCpL, TL)
    );

    // energy equation
    fvScalarMatrix TVEqn
    (
        fvm::ddt(rhoCpV,TV)
        + fvm::div(rhoCpVPhi,TV)
	    - fvm::Sp(fvc::ddt(rhoCpV) + fvc::div(rhoCpVPhi), TV)
        - fvm::laplacian(kappaVEff,TV)
        - dpdt*alpha2
        + phaseChange.TSource2()
        - latentSourceV
            ==
        fvOptions(rhoCpV, TV)
    );

    /*
    label liquidCellSize = 0;
    label gasCellSize = 0;

    labelList FluidStatus(mesh.nCells(),0);

    forAll(alpha1,cellI)
    {
        if (alpha1.oldTime()[cellI] >= (1-cutOff))
        {
            FluidStatus[cellI]= 1;
            liquidCellSize++;
        }
        else if (alpha1.oldTime()[cellI] <= cutOff)
        {
            FluidStatus[cellI]= -1;
            gasCellSize++;
        }
        else
        {
            FluidStatus[cellI]= 0;
        }
    }

    labelList  GasCells(gasCellSize);
    scalarList GasSatTemp(gasCellSize);

    labelList  LiquidCells(liquidCellSize);
    scalarList LiquidSatTemp(liquidCellSize);

    label liquidCellCounter = 0;
    label gasCellCounter = 0;

    forAll(FluidStatus,cellI)
    {
        if (FluidStatus[cellI]== -1) // gas
        {
            GasCells[gasCellCounter] = cellI;
            GasSatTemp[gasCellCounter] = TSat[cellI];
            gasCellCounter++;
        }

       // if (FluidStatus[cellI] == 1) // liquid
       // {
       //     LiquidCells[liquidCellCounter]= cellI;
       //     LiquidSatTemp[liquidCellCounter] = TSat[cellI];
       //     liquidCellCounter++;
       // }
    }
    */

    TL = max(min(TL,TMax),TMin);
    TV = max(min(TV,TMax),TMin);

    //TLEqn.setValues(GasCells, GasSatTemp); // gas cells are set to sat Temp
    //TVEqn.setValues(LiquidCells, LiquidSatTemp); // liquid cells are set to sat temp
   
    TL.correctBoundaryConditions();
    TV.correctBoundaryConditions();

    TLEqn.relax();
    TVEqn.relax();

    TLEqn.solve();
    TVEqn.solve();

    // --- 相间热平衡 "整理" 步骤 ---
    // 计算每个网格的总显热 J
    volScalarField totalSensibleHeat = rhoCpL * TL * alpha1 + rhoCpV * TV * alpha2;
    // 计算每个网格的总热容 J/K
    volScalarField totalHeatCapacity = rhoCpL * alpha1 + rhoCpV * alpha2;
    // 计算平衡温度 K
    volScalarField TEquilibrium = totalSensibleHeat / totalHeatCapacity;

    // loose control
    // TV = ((1-omega_) * rhoCpL * alpha1 * TV + omega_ * rhoCpV * alpha2 * TL )/ ((1-omega_) * rhoCpL * alpha1 + omega_ * rhoCpV * alpha2 );
    TV = TL;

    // 更新边界条件
    TL.correctBoundaryConditions();
    TV.correctBoundaryConditions();

    TL = max(min(TL,TMax),TMin);
    TV = max(min(TV,TMax),TMin);

    Info << "min(TL) " << gMin(TL) << "  max(TL) " << gMax(TL)  << endl;
    Info << "min(TV) " << gMin(TV) << "  max(TV) " << gMax(TV) << endl;

    T = TEquilibrium;

    T.correctBoundaryConditions();

    mixture.correctThermo(TL,TV);

    mixture.correct();

    Info<< "min(T) " << min(T).value() << endl;

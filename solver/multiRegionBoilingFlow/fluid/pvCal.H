 //phaseChange.correctSatProperties(p,T);
    // 用水蒸气分压计算饱和温度
    // p_vapour = Y_vapour * p_total * (M_air/M_vapour)
    // 对于理想气体混合物：p_vapour = x_vapour * p_total（摩尔分数）
    // 简化：用质量分数近似
PtrList<volScalarField>& Y2 = mixture.phase2().Y();

    forAll(Y2, i)
    {
        if (Y2[i].name() == "vapour.gas")  // 注意带phase后缀
        {
            p_vapour = Y2[i] * p * scalar(28.97/18.015);
            break;
        }
    }

//scalar cutOff = 1e-3;

forAll(alpha1, cellI)
{
    scalar a = alpha1[cellI];

    if (a > (1-cutOff))
    {
        // 纯液相：用总压（液态水饱和压力由T决定，这里用总压近似）
        p_vapour[cellI] = p[cellI];
    }
    else if (a < cutOff)
    {
        // 纯气相：用水蒸气分压
        p_vapour[cellI] = max
        (
            Y2[cellI] * p[cellI] * scalar(28.97/18.015),
            scalar(100.0)  // 最小100Pa防止log(0)
        );
    }
    else
    {
        // 界面：加权
        scalar pV_gas = max
        (
            Y2[cellI] * p[cellI] * scalar(28.97/18.015),
            scalar(100.0)
        );
        p_vapour[cellI] = a * p[cellI] + (1.0-a) * pV_gas;
    }
}

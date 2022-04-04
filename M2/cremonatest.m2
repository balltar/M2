cremona = (I, n) -> (
    g := I_*;
    if n == 3 then ideal(g_1*g_2,g_0*g_2,g_0*g_1)
    else if n == 2 then ideal(g_2*g_2,g_0*g_1,g_0*g_2)
    )

HB = I -> (res I).dd_2

dualParam = I -> (
    G := gens I;
    J := diff(transpose vars (ring I),G);
    I' := (det sub'(J,{0}),-1*det sub'(J,{1}),det sub'(J,{2}));
    g := gcd(I'_1, I'_2, I'_3);
    ideal(I'_1//g,I'_2//g,I'_3//g)
    
restart
needsPackage "ReesAlgebra"
R = QQ[x,y]
Ma = matrix{{x,x^5+y^5},{y,x^5+x*y^4},{0,y^5}}
Mb = matrix{{x^2,x^4+y^4},{y^2,x^4+x*y^3},{0,y^4}}
Mc = matrix{{x^2,x^4+y^4},{y^2,x^4+x*y^3},{x*y,y^4}}
Md = matrix{{x^3+y^3,x*y^2},{x^3+x*y^2,x^3+x^2*y+y^3},{x^3+x^2*y+y^3,y^3}}
Me = matrix{{x^3+y^3,x^3+x^2*y+y^3},{x^3+x*y^2,x^3+y^3},{0,y^3}}
Mf = matrix{{x^3+y^3,x^3+y^3},{x^3+x*y^2,x^3+x^2*y+y^3},{0,y^3}}
Mg = matrix{{x^3+y^3,x^3+y^3},{x^3+x*y^2,0},{0,y^3}}

--Assuming homogeneous
bidegrees = M -> (
    I1    = minors(2, M);
    delta = (degree I1_0)_0;
    J     = reesIdeal I1;
    bloo = apply(J_*, i -> degree i);
    apply(bloo, i -> i-{0,i_0*delta})
    )

bidegrees(Ma)
bidegrees(Mb)
bidegrees(Mc)
bidegrees(Md)
bidegrees(Me)
bidegrees(Mf)
bidegrees(Mg)


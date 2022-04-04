--Given a polynomial ring R=k[R_*], give a basis of degree n part
--of the subring k[variables]. Return this basis as a list of
--elements of R.
subBasis = (n,R,variables) -> (
    local i;
    S := coefficientRing(R)[variables];
    apply(flatten entries basis(n,S),i->sub(i,R))
    )

--generate random homog poly of degree d in R
genHPoly = (R,d) -> (
    sum apply(flatten entries basis(d,R), i->random(coefficientRing(R))*i)
    )

--formats generators so you can copy them into desmos!
desTex = idealGen -> (
    concatenate separate("$",concatenate(separate("\\,",tex idealGen)))
    )

--Take dot product of two lists of same length
dot = (ls1, ls2) -> (
    local i;
    l := length ls1-1;
    sum toList apply(0..l,i->(ls1_i)*(ls2_i))
    )

--Given an element of a polynomial ring R=k[rest,variables] give
--the coefficient vector of it viewing R as k[rest][variables]
coefficientsWRT = (r, variables) -> (
    local i;
    R := ring(r);
    rest := R_*-set(variables);
    S := coefficientRing(R)[rest][variables];
    apply(flatten entries (coefficients(sub(r,S)))_1, i->sub(i,R))
    )

--Given two ring elements r and s outputs the maximum i such that
--r^i divides s.
maxPowerDividing = (r,s) -> (
    i := 0;
    while true do (
	if s%(r^i) == 0 then i = i+1 else return i-1;
	)
    )

--Given a polynomial r, gives the maximum power of the specified
--variable appearing in the terms of r 
maxDegInTerms = (r,variable) -> (
    local t;
    max apply(terms r,t -> maxPowerDividing(variable, t))
    )

--Given a polynomial r, outputs which variables appear with maximum
--power of one in the terms of r.
linearVariables = r -> (
    local v;
    select((ring r)_*, v->maxDegInTerms(r,v)==1)
    )

--Given a homoegenous ring element r, outputs list of full coefficients
--with respect to basis (so includes zeroes)
fullCoefficients = r -> (
    d := first degree r;
    R := ring(r);
    basisMons := flatten entries basis(d,R);
    if (zero r) then return apply(toList(1..length basisMons), i->substitute(0,R))
    else return for i in basisMons list coefficient(i,r);
    )

--Given an element r of a polynomial ring k[rest,variables],
--constructs a matrix by viewing r as an element of 
--k[rest][variables], extracting the linear variables from rest,
--and then making equations coefficientwise, viewing r as an element
--of k[rest][variables]
linearExtractMatrix = (r,variables) -> (
    local i;
    R := ring r;
    kk := coefficientRing(R);
    rest := R_*-set variables;
    S  := kk[rest][variables];
    coeffs := flatten entries (coefficients(sub(r,S)))_1;
    --Put something to actually verify no two linear terms in same coefficient
    use R;
    lv := (linearVariables r)-set variables;
    restv := R_* - (set variables + set lv);
    Sv := kk[restv][lv];
    rows := apply(coeffs, i->fullCoefficients(sub(i,Sv)));
    (subBasis(first degree (sub(r,S)),R,variables),lv,mutableMatrix sub(matrix rows,R))
    )

permuteList = (ls, perm) -> (
    apply(perm, i->ls_i)
    )

savedReduce = M -> (
    rowAdd(
	rowAdd(
	    rowAdd(
		rowAdd(
		    rowAdd(
			rowAdd(
			    rowAdd(
				rowAdd(
				    rowAdd(
					rowMult(
					    rowAdd(
						rowAdd(
						    rowMult(M,1,1/2)
						    ,0,-1,1)
						,2,-1,1)
					    ,2,-2)
					,0,1/2,2)
				    ,1,-1/2,2)
				,0,2,3)
			    ,1,-2,3)
			,2,3,3)
		    ,0,-3,4)
		,1,3,4)
	    ,2,-4,4)
	,3,-2,4)
    )

subList = (ls,ols) -> (
    apply(ls, i->sub(i,ols))
    )

--RUN THE FOLLOWING CODE AND FUNCTIONS ABOVE TO GENERATE AN EXAMPLE
R = QQ[a_1..a_6,b_1..b_6,c_1..c_6,alpha_1,alpha_2,beta_1,beta_2,g_7..g_10,d_7..d_10,x,y]
for i from 1 to 6 do Q_i = dot({a_i,b_i,c_i},subBasis(2,R,{x,y}))
use R;
varphi = matrix{{x^2*Q_1,x^2*Q_2},{(x+y)^2*Q_3,(x+y)^2*Q_4},{y^2*Q_5,y^2*Q_6}}
O = join(newopts,opts,opts2,newopts2)
phi = sub(varphi,O)

--BELOW IS A DERIVATION OF O
for i from 1 to 2 do l_i = dot({alpha_i,beta_i},subBasis(1,R,{x,y}))
for i from 7 to 10 do L_i = dot({g_i,d_i},subBasis(1,R,{x,y}))
use R;
p1 = matrix{{1_R,1_R,0_R}}
p2 = matrix{{0_R,1_R,1_R}}
eqn_1 = (flatten entries(p1*varphi))_0-l_1^3*L_7
eqn_2 = (flatten entries(p1*varphi))_1-l_1^3*L_8
eqn_3 = (flatten entries(p2*varphi))_0-l_2^3*L_9
eqn_4 = (flatten entries(p2*varphi))_1-l_2^3*L_10
for i from 1 to 4 do (M_i = linearExtractMatrix(eqn_i,{x,y}); use R;)
--the permutations for eqn_2 and eqn_4 is
use R
M = savedReduce columnSwap((linearExtractMatrix(eqn_1,{x,y}))_2,4,5)
use R
Mheads = permuteList((linearExtractMatrix(eqn_1,{x,y}))_1,{0,1,2,3,5,4,6,7})
colPerm = {5,4,3,2,0,1,7,6}
rowPerm = {4,3,2,1,0}
use R
N = savedReduce rowPermute(columnPermute((linearExtractMatrix(eqn_3,{x,y}))_2,0,colPerm),0,rowPerm)
use R
Nheads = permuteList((linearExtractMatrix(eqn_3,{x,y}))_1,colPerm)

I1 = ideal(dot(drop((entries M)_1,5),drop(apply(Mheads,i->-1*i),5))-dot(drop((entries N)_4,5),drop(apply(Nheads,i->-1*i),5)))
I2 = ideal(dot(drop((entries M)_3,5),drop(apply(Mheads,i->-1*i),5))-dot(drop((entries N)_3,5),drop(apply(Nheads,i->-1*i),5)))
I3 = ideal(dot(drop((entries M)_4,5),drop(apply(Mheads,i->-1*i),5))-dot(drop((entries N)_1,5),drop(apply(Nheads,i->-1*i),5)))

opts = {g_7=>3,d_7=>1,g_9=>1,d_9=>6,alpha_1=>2,alpha_2=>-1,beta_1=>-1,beta_2=>3,c_1=>-2,a_5=>19}

sub(dot(drop((entries M)_0,5),drop(apply(Mheads,i->-1*i),5)),opts)
sub(dot(drop((entries M)_1,5),drop(apply(Mheads,i->-1*i),5)),opts)
sub(dot(drop((entries M)_2,5),drop(apply(Mheads,i->-1*i),5)),opts)
sub(dot(drop((entries M)_3,5),drop(apply(Mheads,i->-1*i),5)),opts)
sub(dot(drop((entries M)_4,5),drop(apply(Mheads,i->-1*i),5)),opts)

sub(dot(drop((entries N)_0,5),drop(apply(Nheads,i->-1*i),5)),opts)
sub(dot(drop((entries N)_1,5),drop(apply(Nheads,i->-1*i),5)),opts)
sub(dot(drop((entries N)_2,5),drop(apply(Nheads,i->-1*i),5)),opts)
sub(dot(drop((entries N)_3,5),drop(apply(Nheads,i->-1*i),5)),opts)
sub(dot(drop((entries N)_4,5),drop(apply(Nheads,i->-1*i),5)),opts)

newopts = {a_1 => 25, a_3 => -1, b_1 => -31, b_3 => 5, c_3 => -1, c_5 => 163, c_3 => -1, b_5 => -138, b_3 => 5, a_3 => -1}
varphi' = sub(varphi, join(opts,newopts))
factor (flatten entries (p1*varphi'))_0
factor (flatten entries (p2*varphi'))_0

opts2 = {g_8 => 1, d_8=>1, g_10 => 1, d_10 => 4, alpha_1=>2, alpha_2=>-1,beta_1=>-1, beta_2=>3,c_2=>-18,a_6=>-3}

use R
M2 = savedReduce columnSwap((linearExtractMatrix(eqn_2,{x,y}))_2,4,5)
use R
Mheads2 = permuteList((linearExtractMatrix(eqn_2,{x,y}))_1,{0,1,2,3,5,4,6,7})
colPerm = {5,4,3,2,0,1,7,6}
rowPerm = {4,3,2,1,0}
use R
N2 = savedReduce rowPermute(columnPermute((linearExtractMatrix(eqn_4,{x,y}))_2,0,colPerm),0,rowPerm)
use R
Nheads2 = permuteList((linearExtractMatrix(eqn_4,{x,y}))_1,colPerm)
sub(dot(drop((entries M2)_0,5),drop(apply(Mheads2,i->-1*i),5)),opts2)
sub(dot(drop((entries M2)_1,5),drop(apply(Mheads2,i->-1*i),5)),opts2)
sub(dot(drop((entries M2)_2,5),drop(apply(Mheads2,i->-1*i),5)),opts2)
sub(dot(drop((entries M2)_3,5),drop(apply(Mheads2,i->-1*i),5)),opts2)
sub(dot(drop((entries M2)_4,5),drop(apply(Mheads2,i->-1*i),5)),opts2)

sub(dot(drop((entries N2)_0,5),drop(apply(Nheads2,i->-1*i),5)),opts2)
sub(dot(drop((entries N2)_1,5),drop(apply(Nheads2,i->-1*i),5)),opts2)
sub(dot(drop((entries N2)_2,5),drop(apply(Nheads2,i->-1*i),5)),opts2)
sub(dot(drop((entries N2)_3,5),drop(apply(Nheads2,i->-1*i),5)),opts2)
sub(dot(drop((entries N2)_4,5),drop(apply(Nheads2,i->-1*i),5)),opts2)

newopts2 = {a_2 => 9, a_4 => -1, b_2 => -9, b_4 => 7, c_4 => -1, c_6=> 109, c_4=>-1, b_6 =>-86,b_4 =>7,a_4 =>-1}
varphi'' = sub(varphi, join(opts2,newopts2))
factor (flatten entries (p1*varphi''))_1
factor (flatten entries (p2*varphi''))_1

O = join(newopts,opts,opts2,newopts2)
phi = sub(varphi,O)

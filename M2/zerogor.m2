--Conventions: Throughout n will denote the number of variables,
--kk will denote the base field
--d will denote the power of maximal ideal/degree an ideal is equigenerated in
--R will be a polynomial ring
--I will be an ideal 
--J will be its ReesIdeal
--f will be a general element of a ring
--ci will be a complete intersection ideal
--A, M, C, T will be matrices
--t will be variable in hilbert series
--x_i will be variables
--m will be an integer such that (x_1^m..x_n^m) is contained in I
--mm will be the homogeneous maximal ideal
--F(I) will be the special fiber ring of I
--S will be the polynomial ring for the numerator of the hilbert series
--i,j,k will be auxillary variables used in loops and apply statements

--create a random homogeneous polynomial of degree d given ring R.
--Has option to have coefficients be integers
randPoly = {int => false} >> opts -> (d,R) -> (
    kk := coefficientRing(R);
    if opts.int then
    sum(apply(flatten entries basis(d,R),i->random(ZZ)*i))
    else sum(apply(flatten entries basis(d,R),i->random(coefficientRing R)*i))
)

--Computes minimal number of generators of the ideal
mu = I -> length flatten entries mingens I

--Computes a random artinian Gorenstein-linear ideal I in n variables
--generated in degree d such that I contains (x_1^m..x_n^m).
--Notice that any artinian gorenstein ideal contains (x_1^m..x_n^m) for some m
gorLin = {special => false} >> opts -> (d,n,m) -> (
    if(m<d) then error "There are no homogeneous artinian Gorenstein-linear ideals with d<m!";
    R := QQ[symbol x_1..symbol x_n];
    e := n*(m-1)-2*(d-1);
    --print e;
    ci := ideal apply(R_*,i->i^m);
    local f;
    if opts.special then f = (sum R_*)^e else f = randPoly(e,R, int=>true);
    ideal(flatten entries mingens(ci:ideal f))
    )

randomGorForm = (d,n,m) -> (
    if(m<d) then error "There are no homogeneous artinian Gorenstein-linear ideals with d<m!" ;
    e := n*(m-1)-2*(d-1);
    R := QQ[symbol x_1..symbol x_n];
    f := randPoly(e,R,int=>true)
    )

--gorLinGenMat makes matrix associated to directrix
--Assumes polynomial ring of f is standard graded
gorLinMat = {testing => false, generic => "false", directrix => 0} >> opts -> (delta,n,m) -> (
    e := n*(m-1)-2*(delta-1);
    S := QQ[symbol x_1..symbol x_n];
    local f,R;
    if zero opts.directrix then (
    if opts.generic === "true" then (
	B := flatten entries basis(e,S/(ideal(S_*))^[m]);
	lB := length B;
	R = frac(QQ[symbol a_0..symbol a_(lB-1)])[symbol x_1..symbol x_n];
	f = sub(sum toList apply(0..lB-1, i-> a_i*sub(B_i,R)),R);
	)
    else (
	if opts.generic === "special" then (
	    R = frac(QQ[symbol a_1..symbol a_n])[symbol x_1..symbol x_n];
	    f = (sum toList apply(1..n, i-> a_i*x_i))^e;
	    )
	else (
	    R = S;
	    f = randPoly(e,R,int=> false);
	);
    );
    )
    else (
	f = sub(opts.directrix, S);
	);
    kk := coefficientRing R;
    gamma := flatten entries sub(basis(e+delta,R/((ideal(R_*))^[m])),R);
    beta := flatten entries basis(delta,R);
    M := matrix toList apply(gamma,i->toList apply(beta,j->(alpha := i//j; if zero alpha then 0_kk else f_(alpha))));
    if opts.testing then (
    	print("e equals",e);
    	print("f equals",f);
    	print("gamma is", gamma);
    	print("beta is", beta);
	K := sub(generators ker M, R);
	I := ideal (basis(delta,R)*K);
    	J := ((ideal(R_*))^[m]):ideal(f);
    	print I;
    	print J;
    	print(I == J);
    	);
    M
    )

completeSymmetricIdeal = (n,m,d) -> (
    R := QQ[symbol x_1..symbol x_n];
    ((ideal R_*)^[m]):ideal sum flatten entries basis(n*(m-1)-2*(d-1),R)
    )

compSym = n -> (
    R := QQ[symbol x_1..symbol x_n];
    ((ideal R_*)^[m]):ideal sum flatten entries basis(n*(2*n-2-1)-2*(n-1),R)
    )

--Newton Duals
newtonMatrix = f -> (
    matrix apply(exponents f, i->vector i)
    )

coeffs = f -> (
    flatten entries (coefficients f)_1
    )

expToMons = (v, R) -> (
    v' := entries v;
    product toList apply(0..length v'-1, i->(sub((R_*)_i,frac R))^(v'_i))
    )

newtonMatrixToPoly = (f, M) -> (
    R := ring f;
    C := coeffs f;
    sum toList apply(0..rank source M-1, i-> C_i*expToMons(M_i, R))
    )

newtonDual = (f,m) -> (
    alpha := vector toList apply(0..length (ring f)_*-1, i-> m-1);
    A := matrix toList apply(0..length coeffs f -1, i->alpha);
    --print A;
    newtonMatrixToPoly(f,A - newtonMatrix(f))
    )


--Pure Power Gap


--Express HS(F(mm^d),t) as Q(t,d)/(1-t)^n. 
--Then this computes Q(t,d).  
hVecMaxPower = n -> (
    R := QQ[symbol x];
    A := mutableMatrix sub(((matrix vector join({1_R},apply(1..n-1, i-> 0_R))) | matrix toList apply(1..n-1, i-> vector join(reverse apply(flatten entries (coefficients(product toList apply(1..i,j->x+j_R))_R)_1, k->k*(1/i!)),toList apply(1..n-1-i,j->0)))),QQ);
    S := QQ[symbol d][symbol t];
    M := mutableMatrix sub(matrix toList apply(0..n-1, i->vector matrix solve(A,(mutableMatrix id_(QQ^n))_{i})),S);
    C := mutableMatrix aFCoeff(binomial(d*t+n-1,n-1));
    T := mutableMatrix matrix {reverse toList apply(0..n-1,i->(1-t)^(i))};
    (T*M*C)_(0,0)
    )

--Express HS(F(mm^d),t) as Q(t,d)/(1-t)^n. 
--Then this computes Q(t,d').  
hvemp = (n,d') -> (
    sub(hVecMaxPower n,{d=>d'})
    )

--Checks bidegrees of rees algebra
bidegreesRees = I -> (
    J := reesIdeal I;
    apply(J_*,i->degree i)
    )

--Used for detecting which orbit in section 7 of first Kustin Khoulray paper
factorGens = I -> apply(I_*, i->((factor i)#0)#1)


--The following is an example of an f in the unusual Kustin class:
  5     4       3 2    2 3       4     5     4       3         2 2         3       4       2   2       2 2     3 2     2 3         3     2 3       4     5
2x  + 2x x  + 2x x  + x x  + 6x x  + 7x  + 8x x  + 8x x x  + 9x x x  + 7x x x  + 2x x  + 3x x x  + 4x x x  + 9x x  + 3x x  + 8x x x  + 8x x  + 6x x  + 4x
  1     1 2     1 2    1 2     1 2     2     1 3     1 2 3     1 2 3     1 2 3     2 3     1 2 3     1 2 3     2 3     1 3     1 2 3     2 3     2 3     3
R = QQ[symbol x_1, symbol x_2, symbol x_3];
mm = ideal R_*;
d = 5;
f = 2*x_1^5+2*x_1^4*x_2+2*x_1^3*x_2^2+x_1^2*x_2^3+6*x_1*x_2^4+7*x_2^5+8*x_1^4*x_3+8*x_1^3*x_2*x_3+9*x_1^2*x_2^3*x_3+7*x_1*x_2^3*x_3+2*x_2^4*x_3+3*x_1^2*x_2*x_3^2+4*x_1*x_2^2*x_3^2+9*x_2^3*x_3^2+3*x_1^2*x_3^3+8*x_1*x_2*x_3^3+8*x_2^2*x_3^3+6*x_2*x_3^4+4*x_3^5;
I = (mm^[d]):ideal(f)

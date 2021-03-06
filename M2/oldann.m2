--Make a linear combination of elements of lsvars with random coefficients from coefficientRing of testRing above
randLinComb = (lsvars) -> (
        sum apply(lsvars,z-> random(coefficientRing(ring(z)))*z)
    )

--Given a matrix and a subset of its entries, replace all entries by a subset ls of entries.
--Notice you can't replace zero.
specializeVars = (R,A,ls) -> (
    nn := numgens target A - 1;
    M := new MutableMatrix from A;
    for i from 0 to nn do (
	for j from i+1 to nn do (
	    if not member(M_(i,j), ls) and M_(i,j)!=0 then (
		r:=randLinComb(ls);
		M_(i,j)=r;
		M_(j,i)=-r;
		)
	    )
	);
    matrix M
    )

--Randomly specialize variables in matrix so all entries are linear combinations of n variables
randSpecializeVars = {showBasis => false} >> opt -> (variables,A,n) -> (
    rls := take(apply(random(for i from 1 to length variables-1 list i),  i->R_i), n);
    if opt.showBasis then (reverse sort rls, specializeVars(R,A,rls)) else specializeVars(R,A,rls)
    )

--Make a generalSkewMatrix
generalSkewMatrix = {showBasis => false} >> opt -> (R,n,m) -> (
    if m>n*(n-1)//2 then error "Not enough variables to make that specialization";
    randSpecializeVars(showBasis => opt.showBasis, for i from 0 to n*(n-1)//2-1 list R_i, genericSkewMatrix(R,R_0,n),m)
    )

--Make a general BMatrix (can specify which matrix in resolution you want
generalBMatrix = {showBasis => false, whichD => 2} >> opt -> (R,n,m) -> (
    if opt.showBasis then (
	output := generalSkewMatrix(R,n,m,showBasis=>true);
	(output_0, genericBMatrix(n,whichD=>opt.whichD,useSkewMatrix=>output_1))
	)
    else genericBMatrix(n,whichD=>opt.whichD,useSkewMatrix=>generalSkewMatrix(R,n,m))
    )

--Make a general BIdeal from ring, dimensions matrix, number of vars to keep
generalBIdeal = {showBasis => false} >> opt -> (R, n, m) -> (
    if opt.showBasis then (
	output = generalBMatrix(R,n,m,whichD => 1,showBasis=>true);
	(output_0, ideal output_1)
	)
    else ideal generalBMatrix(R,n,m,whichD => 1)
    )

monomialsFromIdeal = (coeffRing,J) -> (
	for g in J_* list (
	    C := flatten entries (coefficients g)_1;
	    for c in C list flatten entries monomials substitute(c, coeffRing)
	    )
    )

--outputs true if the monomials are as expected
expectedMonomial = () -> (
    use R;
    I := genericBIdeal(6,useSkewMatrix => specializeVars(R,genericSkewMatrix(R,6),{b,c,d,e}));
    print {b,c,d,e};
    M := monomialsFromIdeal(R, reesIdeal I);
    count := 0;
    for s in subsets(for i from 1 to length R_*-3 list R_i, 4) do (
	use R;
	s={b,c,h,j};
	II := genericBIdeal(6,useSkewMatrix => specializeVars(R,genericSkewMatrix(R,6),s));
	print s;
	count = count +1;
	print count;
	use R;
	MM := monomialsFromIdeal(R, reesIdeal substitute(II, {s_0=>b,s_1=>c,s_2=>d,s_3=>e}));
	print(M == MM);
	print "#####################";
	)
    )

reesTest = {checkLinear => false} >> opt -> (R,n,m) -> (
    	I = generalBIdeal(R,n,m);
	J = reesIdeal I;
        if opt.checkLinear then print concatenate("Is linear type? ",toString((isLinearType I)));
	for i in J_* list degree i
    )

 R = QQ[vars(0..9),z_1,z_2]
 A = genericSkewMatrix(R,R_0,5)
 R = QQ[vars(0..14),z_1,z_2]
 A = genericSkewMatrix(R,R_0,6)


--Generic | Skewsym tests
--ed2g = (n,k,t,R) -> (
--    	A := genericMatrix(R,n,k);
--	B := genericSkewMatrix(R,R_(k*n),n);
--	M := A | B;
--	--(M,concatenate("Codim of ",toString(t),"x",toString(t)," minors is: "),codim minors(t, M))
--	codim minors(t, M)
--    )

--The following illustrates some weird behavior:
--submatrix'(T,{-1},{-1})

--Creates a vector of heights of ith Fitting ideal for i=1,2,...,max size minor-1 to check G_d
gVector = A -> (
    nn := min {numgens target A, numgens source A};
    for i from 1 to nn-1 list codim minors(nn-i,A)
    )

isGInf = A -> (
    nn := min {numgens target A, numgens source A};
    gvec := for i from 1 to nn-1 list (print concatenate("Computing F_",toString i);codim minors(nn-i,A));
    all(nn-1,i-> i<=gvec_i)
    )

resFormat = C -> (for i from 0 to length C list rank C#i)

--Specialize 5x5 to 5 vars is linear type






--HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRREEEEEEEE
--HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRREEEEEEEE
--HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRREEEEEEEE
--HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRREEEEEEEE
--HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRREEEEEEEE
--HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRREEEEEEEE

genericBMatrix = {useRing => ZZ/9973, useSkewMatrix => {}, usez1 => {}, usez2 => {}, whichD => 2, reduceEntries => false} >> opt -> n -> (
    R  := opt.useRing;
    A  := opt.useSkewMatrix;
    z1 := opt.usez1;
    z2 := opt.usez2;
    if R === ZZ/9973 then R = testRing(n);
    if A === {} then A = genericSkewMatrix(R,n) else R = ring(A_(0,0));
    if numgens target A != n then error "Dimension of matrix isn't same as n";
    if z1 == {} then z1 = R_(n*(n-1)//2);
    if z2 == {} then z2 = R_(n*(n-1)//2+1);
    if opt.reduceEntries then (
	M = new MutableMatrix from A; M_(0,1)=0_R; M_(1,0)=0_R; A = matrix M;
	); 
    print A;
    nn   := n - 1;
    --First case is even
    if nn%2 == 1 then (
	else if opt.whichD == 2 then (
    	     B1  := matrix {{-subpf(A,{0,1}), z1, z2}, {pf(A),z1*A_(0,1),z2*A_(0,1)}};
    	     B2  := matrix {for j from 2 to nn list 0, for j from 2 to nn list z2*A_(0,j)-z1*A_(1,j)};
    	     (B1 | B2) || ((matrix for i to nn-2 list {0}) | transpose A_{2..nn})
	    )
	else if opt.whichD == 3 then (
    	    ls1  := {-z1,-subpf(A,{0,1}),0_R}|for j from 2 to nn list subpf(A, {1,j});
    	    ls2  := {z2,0,subpf(A,{0,1})}|for j from 2 to nn list subpf(A, {0,j});
    	    transpose matrix {ls1,ls2}
	    )
	else error "expected 1, 2, or 3 to specify which differential"
	)
    --Next case is odd
    else (
	else if opt.whichD == 2 then (
	    toprow := matrix {join({-subpf(A,{1}), subpf(A,{0})}, for i to nn-2 list 0)};
    	    A' := matrix join({{0, z2*A_(0,1)-z1}, {-z2*A_(0,1)+z1,0}}, for i to nn-2 list {-z2*A_(0,i+2), -z2*A_(1,i+2)});
    	    transpose (toprow || (A' | submatrix(A,{2..nn})))
    	    )
	else if opt.whichD == 3 then (
	    topblock    := matrix{{z1,-z2},{subpf(A,{0}),0}, {subpf(A,{1}),0}}; 
    	    bottomblock := matrix for i from 2 to nn list {subpf(A,{i}),subpf(A,{0,1,i})};
    	    topblock || bottomblock
	    )
	else error "expected 1, 2, or 3 to specify which differential"
	)
    )

genericBComplex = {useRing => ZZ/9973, useSkewMatrix => {}} >> opt -> n -> (
	R := opt.useRing;
    	if R === ZZ/9973 then R = testRing(n);
	C      := new ChainComplex;
	C.ring  = R;
	M1     := genericBMatrix(useRing => R,useSkewMatrix => opt.useSkewMatrix, n,whichD => 1);
	M2     := genericBMatrix(useRing => R,useSkewMatrix => opt.useSkewMatrix, n,whichD => 2);
	M3     := genericBMatrix(useRing => R,useSkewMatrix => opt.useSkewMatrix, n,whichD => 3);
	C#0    = target M1;
	C#1    = source M1;
	C#2    = source M2;
	C#3    = source M3;
	d1     := map(C#0,C#1,M1);
	d2     := map(C#1,C#2,M2);
	d3     := map(C#2,C#3,M3);
	C.dd#1 = d1;
	C.dd#2 = d2;
	C.dd#3 = d3;
	if d1*d2 == 0 then (if d2*d3 == 0 then C else error "d_2*d_3 is not zero")
	else error "d_1*d_2 is not zero"
    )

--Specialize variables in generators of an ideal
specializeIdeal = {rand => false} >> opt -> (I,n) -> (
    local vls;
    R := ring I;
    if opt.rand then vls = take(random R_*, n) else vls = apply(n,i-> R_i);
    T := (coefficientRing R)[vls];
    M := map(T,R,matrix{join(T_*,apply(length R_* - n, i -> randLinComb(T)))});
    M(I)
    )
restart

--pfaffian of matrix
pf = A -> (pfaffians(numgens target A, A))_0

--Compute the sign in the signed pfaffian
pfsgn = ls -> ((-1)^(sum(ls)%2+sum(for j from 0 to length ls-1 list sum(for i from 0 to j list if ls_i>ls_j then 1 else 0))))

--take a pfaffian of submatrix by deleting rows and columns in index list ls
subpf = (A,removedIndices)     -> (
    nn := numgens target A - length removedIndices;
    if odd nn then 0
    else (
	if nn == 0 then 1
	else pfsgn(removedIndices)*(pfaffians(nn, submatrix'(A,removedIndices,removedIndices)))_0)
    )

--n*(n+1)//2 is the nth triangle number, i.e., number of vars in skewsym matrix
testRing = n -> (QQ[symbol x_1..symbol x_((n*(n+1))//2)])

--Makes a random linear combination of variables
randLinComb = ls -> sum apply(length ls,i->sub(random coefficientRing (ring (ls_0)),R)*ls_i)

generateGeneric = (n,m) -> (
    R := QQ[symbol t_1..symbol t_((n*(n+1))//2),symbol z_1,symbol z_2];
    T := genericSkewMatrix(R,n);
    shuffle := random(toList(1..(n*(n+1))//2));
    rvars := apply(toList(1..m), i-> (R_*)_(shuffle_i))
    varis := set R_* - rvars;
    T' = sub(T,apply(rvars, i-> i => randLinComb varis))
    use R;
    ideal(pf(T), subpf(T,{0,1}), apply(2..(n-1),i->sub(randLinComb(R'),R)*subpf(T,{0,i})+sub(randLinComb(R'),R)*subpf(T,{1,i})))
    )

generateGeneral = n -> (
    R := QQ[symbol t_1..symbol t_((n*(n+1))//2),symbol z_1,symbol z_2];
    T := genericSkewMatrix(R,n);
    R' := QQ[delete(R_*,{z_1,z_2})]
    T = sub(T,{z_1 => randLinComb}
    print T;
    ideal(pf(T), subpf(T,{0,1}), apply(2..(n-1),i->z_1*subpf(T,{0,i})+z_2*subpf(T,{1,i})))
    )

--Number of generators of I
mu = I -> (length flatten entries mingens I)

--Checks to see if element of poly ring is a scalar
isScalar = f -> liftable(f, coefficientRing ring f)

--Specialize one at a time, assumes R = k[x_1,...,x_n]. Assumes n<m
specializeMatrix = {rand => false,showMat=>false} >> opt -> (M,n) -> (
    local i,j,X;
    R  := ring source M;
    if opt.rand then X = reverse sort take(random(R_*),n) else X = apply(n,i->R_i);
    T  := (coefficientRing R)[X];
    M' := new MutableMatrix from M;
    N  := new MutableMatrix from substitute(M,T);
    for i to numrows M'-1 do (
	for j from i+1 to numcols M'-1 do (
	    if not member(sub(M'_(i,j),T), T_*) and not isScalar M'_(i,j) then (
		r:=randLinComb(T);
		N_(i,j)=r; N_(j,i)=-r;)));
    if opt.showMat then (print "The specialized matrix is: "; print N;);
    matrix N
    )

--n is dimensions of skew matrix, m is number of variables, R is coefficient ring 
generalSkewMatrix = (R,n,m) -> (
    T := R[x_1..x_m];
    M := mutableIdentity(R,n);
    for i from 0 to n do (
	for j from i to n do (
	    if i=j then M_(i,j)=0
	    else (
		r := randLinComb(T);
		M_(i,j)=r; M_(j,i)=-r;
		);
	    );
	);
    matrix M
    )

--Change so z_1 and z_2 are general linear combs (if only change entries of matrix to be 
--general lin combs of other entries of matrix and don't specialize z_i then not G_d)
--Generic Anne Brown Ideal
genericBIdeal = {useMatrix => matrix {{}}} >> opt -> n -> (
    local M,i;
    if opt.useMatrix == matrix {{}} then M = genericSkewMatrix(testRing(n),n) else M=opt.useMatrix;
    R := ring source M;
    T := (coefficientRing(R))[R_*,symbol z_1,symbol z_2];
    M = substitute(M,T);
    use T;
    if (numgens source M)%2==0 then ideal join({pf(M),subpf(M,{0,1})},for i from 2 to m-1 list z_1*subpf(M,{0,i})+z_2*subpf(M,{1,i}))
    else ideal join({subpf(M,{0}),subpf(M,{1})},for i from 2 to m-1 list z_1*subpf(M,{0,1,i})+z_2*subpf(M,{i}))
    )

genericACI = {useMatrix => matrix {{}}} >> opt -> n -> (
    local M;
    if opt.useMatrix == matrix {{}} then M = genericSkewMatrix(testRing(n),n) else M=opt.useMatrix;
    R := ring source M;
    T := (coefficientRing(R))[R_*];
    M = substitute(M,T); 
    use T;
    if (numgens source M)%2==0 then ideal(pf(M),subpf(M,{0,1}),subpf(M,{0,2}),subpf(M,{1,2}))
    else ideal(subpf(M,{0}),subpf(M,{1}),subpf(M,{2}),subpf(M,{0,1,2}))
    )

--Creates a nxn skew sym matrix with m variables as entries and others lin combs of these
specializedMatrix = (n,m) -> specializeMatrix(rand=>true,showMat=>true,genericSkewMatrix(testRing(n),n),m)

--Creates a B ideal from a specialized matrix
specializedBIdeal = (n,m) -> genericBIdeal(n,useMatrix=>specializedMatrix(n,m))
specializedACI = (n,m) -> genericACI(n,useMatrix => specializedMatrix(n,m))

----------------------------------EXAMPLES------------------------------
--5x5 Case
I = specializedBIdeal(5,12);--Don't specialize any of the 10 variables in the matrix
I = specializedBIdeal(5,10);--Still don't specialize any of the 10 variables in the matrix
I = specializedBIdeal(5,5); --Specialize the matrix to 5 variables, so there are 7 vars with z's
I = specializedBIdeal(5,3); --First interesting case: specialize matrix to 3 vars and z's
--Notice that the nonlinear defining equations only include z_2
isLinearType I --For some reason this runs incredibly slow...
reesIdeal I --Can just see the nonlinear equations

--6x6 Case
I = specializedBIdeal(6,4);
J = reesIdeal I;
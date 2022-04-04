-- Exercise Session 5

Q = QQ[x,y,z]

M = matrix{{x,y,z,0}, {0,x,y,z}}

I = minors(2,M)

betti res I

hvec = I -> (
    ls   = {};
    i    = 0;
    next = hilbertFunction(i,ring(I)/I);
    while not (next == 0) do (
	ls = append(next, ls);
	i  = i+1;
	next = hilbertFunction(i,ring(I)/I);
	);
    ls
    )

hvecN = (I,n) -> for i from 0 to n list hilbertFunction(i,ring(I)/I)


hvec(I,5)

s1=3
s2=3
F1=random(s1,Q)
F2=random(s2,Q)

I=ideal fromDual matrix{{F1,F2}}
 
betti res I
res I

hvec(I,3)

needsPackage "TorAlgebra"
torAlgClass (Q/I)

--I is a homogeneous ideal in a polynomial ring Q such that Q/I is artinian
socdeg = (I) -> (
    	i = 0;
    	while not (hilbertFunction(i,ring(I)/I) == 0) do i=i+1;
	i	
    )

hvecCI = I -> (
    d1 = (degree(I_0))_0;
    d2 = (degree(I_1))_1;
    d3 = (degree(I_2))_2;
    s = d1+d2+d3-3;
    hvec(I,s)
    )

--This can be improved so don't have to compute ideal... just do it directly since CI
hvecLink = (I,d1,d2,d3) -> (
    hvecCI(ideal(x^d1,y^d2,z^d3)) - reverse hvecN(I,d1+d2+d3-3);
    )

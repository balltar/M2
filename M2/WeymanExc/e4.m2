restart
needsPackage "TorAlgebra"
needsPackage "DGAlgebras"

Q = QQ[x,y,z]
I = ideal(x^2,y^2,z^2,y*z)
res I

R = Q/I
A = HH(koszul vars R);
mingens A_1
mingens A_2
mingens A_3

use Q

--notice in following three examples how we are decreasing the betti numbers through links
X = ideal(x^3,y^3,z^3)
J = (X:I)
res J

X = ideal(x^3,y^3,z^2)
J = (X:I)
res J

X = ideal(x^3,y^2,z^2)
J = (X:I)
res J

X = ideal(x^3+y*z,y^2,z^2) --perterb third minimal generator
J = (X:I)
res J
(res J).dd --the culprit here is there's a one in the differential

needsPackage "PruneComplex" 
Qm = localRing(Q, ideal vars Q)
pruneComplex ((res J)**Qm)

use Q
X = ideal(x^2,y^2,z^2)
J = (X:I)
res J

--Exercise 4.2

needsPackage "RandomIdeals"

--Input: generator in a *-local ring and the ring itself
--Output: some perturbation of the ideal using a random monomial ideal
perturbIdeal = (x,Q) -> (
    	x+randomMonomialIdeal({random(10),random(10),random(10)},Q)_0   
    )

--Input: homogeneous ideal
--Output: maximal 
maxreg = (I,Q) -> (
    	g     = codim I;
	isRegular = false;
	while(not isRegular) do (
	    proposedgens = ideal(for i from 0 to g-1 list perturbIdeal(x,Q));
	    if codim proposedgens == g then isRegular = true;
	);
    	proposedgens
    )
    
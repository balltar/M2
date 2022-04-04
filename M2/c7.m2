restart
S = QQ[x,y,z]
T = QQ[t]
K = ker(map(T,S,{t^3,t^4,t^5}))
R = S/K
mm = ideal(x,y,z)
needsPackage "ReesAlgebra"
multiplicity mm

multiplicity3tuple = (a1,a2,a3) -> (
    	Source = QQ[x,y,z];
	Target = QQ[t];
	Kernel = ker(map(Target,Source,{t^a1,t^a2,t^a3}));
	QuotientR = Source/Kernel;
	maxIdeal = ideal(x,y,z);
	multiplicity maxIdeal
    )

----
restart
needsPackage "ReesAlgebra"
computeE = n -> (
    V = QQ[vars(0..2*n-1)];
    X = genericMatrix(V,V_0,2,n);
    I = minors(2,X);
    R = V/I;
    multiplicity ideal(gens R)
    )

for n from 2 to 10 list
    computeE(n)

----
restart
needsPackage "ReesAlgebra"
computeE = n -> (
    V = QQ[vars(0..n-2)];
    Y = matrix{flatten({V_*,0}),flatten({0,V_*})}; 
    I = minors(2,Y);
    R = V/I;
    multiplicity ideal(gens R)
    )

--square of max ideal
for n from 2 to 10 list
    computeE(n)

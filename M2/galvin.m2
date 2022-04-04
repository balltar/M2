restart
R  = QQ[x,y]/(ideal(x*x,y*y))
e1 = x*y
E1 = ideal(e1)
M1 = comodule E1
for i from 1 to 10 list hilbertFunction(i, M1)

----
restart
R  = QQ[x,y,z]/(ideal(x*x,y*y,z*z))
e1 = x*y
e2 = y*z
E1 = ideal(e1)
E2 = ideal(e2)
M1 = comodule E1
M2 = comodule (E1+E2)
for i from 0 to 10 list hilbertFunction(i, M1)
for i from 0 to 10 list hilbertFunction(i, M2)

---
restart
R = QQ[a,b,c,d]
squares = ideal(a*a,b*b,c*c,d*d)
M1 = comodule squares
M2 = comodule (squares + ideal(a*b))
K1 = ker map(M2,M1,matrix{{1_R}})

for i from 0 to 10 list hilbertFunction(i,K1)
for i from 0 to 10 list hilbertFunction(i,M1)
for i from 0 to 10 list hilbertFunction(i,M2)

---
restart
R = QQ[a,b,c,d,e]
squares = ideal(apply(R_*, i-> i^2))
M1 = comodule squares
res M1

--
R  = QQ[a,b,c]
I1 = ideal(a*a,b*b,a*b)
I2 = I1 + ideal(c*c,b*c)
M1 = comodule I1
M2 = comodule I2
K1 = ker map(M2,M1,matrix{{1_R}}) 

S  = QQ[a,b]
J  = ideal(a*a,b*b,a*b)
N  = comodule J
T  = module (S[x]) 
NT = tensor(N,T)

for i from 0 to 10 list hilbertFunction(i,K1)
for i from 0 to 10 list hilbertFunction(i,M1)
for i from 0 to 10 list hilbertFunction(i,M2)

fhilbertFunction = (n, M) -> (
    for i from 0 to n list hilbertFunction(i,M)
    )

--Four generators
restart
needsPackage("GenericInitialIdeal");
R  = QQ[a,b,c,d];
I1 = ideal(a*a, b*b, c*c, a*b, b*c);
I2 = I1 + ideal(d*d, c*d); --path
I3 = I1 + ideal(d*d, b*d); --star

G2 = gin I2
G3 = gin I3
fhilbertFunction(10,comodule I2)
fhilbertFunction(10,comodule I3)
fhilbertFunction(10,comodule G2)
fhilbertFunction(10,comodule G3)

--Five generators
R = QQ[a,b,c,d,e]
four1 = ideal(a*a, b*b, c*c, a*b, b*c)+ideal(d*d,c*d)--path
four2 = ideal(a*a, b*b, c*c, a*b, b*c)+ideal(d*d,b*d)--star
I1 = four1 + ideal(e*e, d*e)--path
I2 = four1 + ideal(e*e, c*e)--fork
I3 = four2 + ideal(e*e, b*e)--star

G1 = gin I1
G2 = gin I2
G3 = gin I3

fhilbertFunction(10, comodule G1)
fhilbertFunction(10, comodule G2)
fhilbertFunction(10, comodule G3)

--Six generators
R = QQ[a,b,c,d,e,f]
five1 = ideal(a*a, b*b, c*c, a*b, b*c)+ideal(d*d,c*d)+ideal(e*e, d*e)--path
five2 = ideal(a*a, b*b, c*c, a*b, b*c)+ideal(d*d,c*d)+ideal(e*e, c*e)--fork
five3 = ideal(a*a, b*b, c*c, a*b, b*c)+ideal(d*d,b*d)+ideal(e*e, b*e)--star

I1 = five1 + ideal(f*f,e*f)
I2 = five1 + ideal(f*f, d*f)
I3 = five1 + ideal(f*f, f*c)
I4 = five2 + ideal(f*f, c*f)
I5 = five2 + ideal(f*f, b*f)
I6 = five3 + ideal(f*f, f*b)

G1 = gin I1
G2 = gin I2
G3 = gin I3
G4 = gin I4
G5 = gin I5
G6 = gin I6

fhilbertFunction(10, comodule G1)
fhilbertFunction(10, comodule G2)
fhilbertFunction(10, comodule G3)
fhilbertFunction(10, comodule G4)
fhilbertFunction(10, comodule G5)
fhilbertFunction(10, comodule G6)

--
R= QQ[a,b,c,d,e,f,g]
I = ideal(a*a,b*b,c*c,d*d,e*e,f*f,g*g)+ideal(a*b,b*c,c*d,d*e,e*f,f*g,g*a)
fhilbertFunction(10, comodule I)

cyclehilbert = (n) -> (
    n' = n-1;
    R = QQ[vars (0..n')];
    I = ideal(for i from 0 to n' list R_i*R_i)+ideal(for i from 0 to n'-1 list R_i*R_(i+1)) + ideal(R_0*R_n');
    print I;
    fhilbertFunction(10,comodule I)
    )

--
restart
fhilbertFunction = (n, M) -> (
    for i from 0 to n list hilbertFunction(i,M)
    )

R = QQ[a,b]
I = ideal(a*a,b*b)
r = a*b
M = comodule I
K = ker(map(M, M, matrix{{r}}))
fhilbertFunction(10,K)
fhilbertFunction(10, comodule I)

J = I + ideal(r)
fhilbertFunction(10,comodule J)

---
isUnimodal = (testls) -> (
    l, i, goinguphill = length testls - 1, 0, true;
    while i < l do (
	if up then (
	    if (testls_i <= testls_(i+1)) then i=i+1
	    else goinguphill = false
	    )
	else if (testls_i >= testls_(i+1) and not goinguphill) then i=i+1
	else return false;
	);
    true
    )

isuhf = (nn,II) -> (
    M = comodule II;
    isUnimodal(for i from 0 to nn list hilbertFunction(i,M))
    )

R = QQ[a,b,c]
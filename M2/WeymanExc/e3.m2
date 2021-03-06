-- Exercise Session 3

needsPackage "TorAlgebra"

Q = QQ[x,y,z]

X = matrix{
    {0,0,0,x,y},
    {0,0,x,y,z},
    {0,-x,0,z,0},
    {-x,-y,z,0,0},
    {-y,-z,0,0,0}}

I = pfaffians(4,X)

res I

for i from 0 to 5 list hilbertFunction(i,Q/I)

isGorenstein (Q/I)

Q = QQ[x,y]

ideal fromDual matrix{{x^2+x*y+y^2}}

F=x^2+x*y+y^2

I = ideal  fromDual matrix{{F}}
It = inverseSystem(matrix{{F}}, DividedPowers => true)
If = inverseSystem(matrix{{F}}, DividedPowers => false)
I == It
I == If

Q=QQ[x,y,z]

I= ideal fromDual matrix{{x^4+x*y*z^2, x*z*y^2}}

res I

s1=4
s2=4
F1=random(s1,Q)
F2=random(s2,Q)

I=ideal fromDual matrix{{F1,F2}}
 
res I

torAlgClass (Q/I)
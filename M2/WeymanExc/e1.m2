--Exercise Session 1--
Q = QQ[x,y,z]

I = ideal(x^2,y^2,z^2,y*z)

F = res I

D = F.dd

r_3 = rank F_3; r_2 = rank F_2 - r_3; r_1 = rank F_1 - r_2;

for i from 1 to 3 list codim minors(r_i,D_i)

K = koszul matrix{{x,y,z}}

A = HH(K)

for i from 0 to 3 list mingens A_i

I = ideal(x^2,y^2,x*z)

codim I

K = koszul matrix{{x^2,y^2,x*z}}

A = HH(K)

for i from 0 to 3 list mingens A_i


Q=QQ[x,y]

I=ideal(x^2*y+y^3, x^3*y, x^3-y^3)

codim I

F=res I

F.dd

Q=QQ[x,y,z]

I=ideal(x^2,x*y,x*z)

codim I

I=ideal(x^2,y^2,x*y)

F=res I

F.dd

Q = QQ[a..h]

--X = genericMatrix(Q,a,2,4)

X=matrix{{a,b,c,d},{e,f,g,h}}

F = eagonNorthcott X

F.dd

Q=QQ[x,y,z]

X = matrix{
    {x,y,z,0},
    {0,x,y,z} }

codim(minors(2,X))

F = res minors(2,X)

F.dd

F = eagonNorthcott X

F.dd
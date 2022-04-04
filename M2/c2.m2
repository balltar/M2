S = QQ[x,y,z,w]
T = QQ[s,t]
I = ker map(T,S,{s^4,s^3*t,s*t^3,t^4})
J = intersect(ideal(x,y),ideal(z,w))
C = S/I
dim C
use S
ideal(x,y,z,w):J
J
a = intersect(I,J)
a:I
--Not CM since linkage preserves CMness and intersection (ideal J above) is not CM)

restart
K = QQ[x,y,z,w]
M1 = matrix{{x,y},{z,w},{x+y,z+w}}
I  = minors(2, M1)
codim I
A = ideal(I_0*x, I_1*y)
J = A:I
A==intersect(J,I)
needsPackage "Depth"
isRegularSequence A_*
I+J
F=res (I+J)
F.dd
codim (I+J)

--
restart
R = QQ[vars (0..7)]
X = genericMatrix(R,R_0,2,4)
I = minors(2,X)
K = koszul(gens I)
H = HH(K)
prune(H)
for i from 0 to 5 list depth(prune HH_i(K))

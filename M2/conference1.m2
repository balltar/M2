R = QQ[x,y,z]/(y^2-x*z)
m = ideal(x,y,z)
S=R[a,b,c,d,e,f]
X = genericMatrix(S,2,3)
a1 = sum(3,i ->  S_i*R_i)
a1 = a*x+b*y+c*z
a2 = d*x+e*y+f*z
J=ideal(a1,a2):m*S
ps = primaryDecomposition J 
apply(ps, x -> codim x)
--compute depth with Ext




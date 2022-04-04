restart
needsPackage "TorAlgebra"

Q = QQ[x,y,z]
I = ideal(x^2,y^2,z^2,y*z)

torAlgClass(Q/I)
torAlgData(Q/I)
torAlgDataList( Q/I , {"m","n","PoincareSeries"} )
torAlgDataPrint( Q/I , {"m","n","PoincareSeries"} )
isCI(Q/I)
isGorenstein(Q/I)
isGolod(Q/I)

P = QQ[a..h]
X = genericMatrix(P,a,2,4)
J = minors(2,X)
torAlgData (P/J)
isGolod (P/J)

needsPackage "DGAlgebras"
torAlgClass (Q/I)
K = koszulComplexDGA (Q/I)
A = homologyAlgebra K

basis(0,A)
basis(1,A)

--Notice in degree two, four of the basis elements are 
--products of basis elements in degree one... i.e.,
--it shows you the small part that has the algebra structure
basis(2,A)
basis(3,A)

Q = QQ[x,y]
torAlgClass (Q/ideal(x^2,y^2))
torAlgClass (Q/ideal(x^2,x*y))

Q = QQ[w,x,y,z]
torAlgClass (Q/ideal(w^2,x^2,y^2,z^2))
I = ideal (x^2,y^2,z^2,x*w,y*w,z*w,w^3-x*y*z)
isGorenstein (Q/I)
torAlgClass (Q/I)
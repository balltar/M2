--Exercise 6.1
--Calculate Buchsbaum-Eisenbud Multipliers
--Calculate Multiplicative Structure
R = QQ[a,b,c,d,e,f,g,h]
M = genericMatrix(R,2,4)
F = eagonNorthcott(M)
D2 = F.dd_2
r2 = rank(F.dd_2)
I2 = minors(r2,D2)

--Buchsbaum-Eisenbud multipliers found by taking minors and then looking for common
--factors across the rows
factor(I2_0)
factor(I2_1)
factor(I2_2)

--To find the equivariant form


D3 = F.dd_3
--Exercise 6.2
--Calculate Buchsbaum-Eisenbud Multipliers
--Calculate Multiplicative Structure
R = QQ[vars (0..8)]
M = genericMatrix(R,3,3)
I = minors(2,M)
F = res I


--Exercise 6.3 
--Calculate Buchsbaum-Eisenbud Multipliers
R = QQ[w,x,y,z]
K = koszul vars R

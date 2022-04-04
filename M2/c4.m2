--part a
V = QQ[vars (0..23)] --Just a bunch of vars
Y = genericMatrix(V,V_0, 3, 3) --3x3 generic 
X = genericMatrix(V,V_9, 3, 1) --3x1 generic
--k[X,Y]
I = minors(1,Y*X)+ideal(det Y)
--Now compute the grade of I:
codim I --Got 3
--Compute how many gens:
length I_* --Got 4
--Need: (alpha_1, alpha_2, alpha_3)=A(I_0,I_1,I_2)
a = minors(1,genericMatrix(V,V_12,3,4)*transpose(gens I))
time L1=a:I
d = length L1_* - codim L1
--

-----------------------
restart
R   = QQ[x1,x2,y1,y2]
Ils = for i from 1 to 6 list (ideal(x1,x2))^i+(ideal(y1,y2))^i
Fls = for i from 0 to 5 list res Ils_i
--Ask about theorem from lecture
g = codim Ils_0
-----------------------
restart
R = QQ[x0,x1,x2,x3]
m = ideal(x0,x1,x2,x3)
needsPackage "RandomIdeals"
Inls = for n from 4 to 7 list
intersect for i from 1 to n list ideal(
    for j from 1 to 3 list (randomElementsFromIdeal({1},m))_0
    )


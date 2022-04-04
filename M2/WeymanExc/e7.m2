restart
needsPackage "TorAlgebra"

-- Type 2

Q = QQ[x,y,z]

a = {y^2};
A = matrix{ 
    {  0  , a_0 }, 
    {-a_0 ,  0  } }
B = matrix{ 
      {z,0}, 
      {x,z},
      {0,-x} }
T = (matrix{{0,0,0},{0,0,0},{0,0,0}}||-transpose B)|(B||A)

I = pfaffians(4,submatrix'(T,{0},{0})) + 
    pfaffians(4,submatrix'(T,{1},{1})) + 
    pfaffians(4,submatrix'(T,{2},{2})) + 	
    pfaffians(2,submatrix'(T,{0,1,2},{0,1,2}))  	

mingens I

torAlgClass(Q/I)

-- Type 3

Q = QQ[x,y,z]

a = {x^2,y^2,z^2}
A = matrix{ 
    {0,a_0,a_1}, 
    {-a_0,0,a_2}, 
    {-a_1,-a_2,0} }
B = matrix{ 
      {x,y,0}, 
      {z,0,x},
      {0,y,z} }

T = (matrix{{0,0,0},{0,0,0},{0,0,0}}||-transpose B)|(B||A)

I = pfaffians(4,submatrix'(T,{0,1},{0,1})) + 
    pfaffians(4,submatrix'(T,{0,2},{0,2})) + 
    pfaffians(4,submatrix'(T,{1,2},{1,2})) + 	
    pfaffians(6,T)

res I
torAlgClass(Q/I)

-- Type 4

Q = QQ[x,y,z]

a = for i from 0 to 5 list random(1,Q);
A = matrix{ 
    {  0  , a_0 , a_1 , a_2}, 
    {-a_0 ,  0  , a_3 , a_4}, 
    {-a_1 ,-a_3 ,  0  , a_5},
    {-a_2 ,-a_4 ,-a_5 ,  0 } }
B = matrix{ 
      {x*y,0,y*z,0}, 
      {0,z*x,0,x*y},
      {0,y*z,x*z,0} }
T = (matrix{{0,0,0},{0,0,0},{0,0,0}}||-transpose B)|(B||A)

I = pfaffians(6,submatrix'(T,{0},{0})) + 
    pfaffians(6,submatrix'(T,{1},{1})) + 
    pfaffians(6,submatrix'(T,{2},{2})) + 	
    pfaffians(4,submatrix'(T,{0,1,2},{0,1,2}));  	

res I

torAlgClass (Q/I)


needsPackage "TorAlgebra"

-- Type 2

Q = QQ[a_1,b_1..b_6]

A = genericSkewMatrix(Q,a_1,2)
B = genericMatrix(Q,b_1,3,2)
T = (matrix{{0,0,0},{0,0,0},{0,0,0}}||-transpose B)|(B||A)

I = pfaffians(4,submatrix'(T,{0},{0})) + 
    pfaffians(4,submatrix'(T,{1},{1})) + 
    pfaffians(4,submatrix'(T,{2},{2})) + 	
    pfaffians(2,submatrix'(T,{0,1,2},{0,1,2}))  	

res I

torAlgClass(Q/I)

-- Type 3

Q = QQ[a_1..a_3,b_1..b_9]

A = genericSkewMatrix(Q,a_1,3)
B = genericMatrix(Q,b_1,3,3)
T = (matrix{{0,0,0},{0,0,0},{0,0,0}}||-transpose B)|(B||A)

I = pfaffians(4,submatrix'(T,{0,1},{0,1})) + 
    pfaffians(4,submatrix'(T,{0,2},{0,2})) + 
    pfaffians(4,submatrix'(T,{1,2},{1,2})) + 	
    pfaffians(6,T)

res I

time torAlgClass (Q/I)


-- Type 4

Q = QQ[a_1..a_6,b_1..b_12]

A = genericSkewMatrix(Q,a_1,4)
B = genericMatrix(Q,b_1,3,4)
T = (matrix{{0,0,0},{0,0,0},{0,0,0}}||-transpose B)|(B||A)

I = pfaffians(6,submatrix'(T,{0},{0})) + 
    pfaffians(6,submatrix'(T,{1},{1})) + 
    pfaffians(6,submatrix'(T,{2},{2})) + 	
    pfaffians(4,submatrix'(T,{0,1,2},{0,1,2}))  	

res I

time torAlgClass (Q/I)

--edge ideal of pentagon is gor grade three given by pfaffians

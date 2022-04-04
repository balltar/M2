restart
R  = QQ[a_1..a_11][x,y];
L3 = a_2*x+a_3*y;
L4 = a_4*x+a_5*y;
L5 = a_6*x+a_7*y;
L6 = a_8*x+a_9*y;
L7 = a_10*x+a_11*y;
phi = matrix {{a_1*x^2*y,x^3},{y^2*L3,y^2*L4},{L7^2*L5,L7^2*L6}}
createSystem(5,3,3,phi)

R = QQ[a_1,a_2,b_1..b_5,c_1..c_5][x,y];
A = a_1*x*y^3+a_2*y^4;
B = b_1*x^4+b_2*x^3*y+b_3*x^2*y^2+b_4*x*y^3+b_5*y^4;
C = c_1*x^4+c_2*x^3*y+c_3*x^2*y^2+c_4*x*y^3+c_5*y^4;
phi = matrix {{x^2,A},{x*y,B},{y^2,C}}
createSystem(5,2,4,phi)

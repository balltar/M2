R = QQ[a,b,c,d,e,m_1..m_4, s_1..s_4, l_1,l_2,t_1,t_2][x,y]
mu = m_1*x^3+m_2*x^2*y+m_3*x*y^2+m_4*y^3
sigma = s_1*x^3+s_2*x^2*y+s_3*x*y^2+s_4*y^3
lambda = l_1*x+l_2*y
tau = t_1*x+t_2*y
phi1 = matrix {{x^2},{y^2},{(x+y)^2}}
phi2 = matrix {{x^2*(x^2+a*x*y+b*y^2)},{y^2*(c*x^2+d*x*y+e*y^2)},{0}}
hxEqns = (diff(x,mu)-sigma)*phi1+(diff(x,lambda)-tau)*phi2+mu*diff(x,phi1)+lambda*diff(x,phi2)
hyEqns = (diff(y,mu)-sigma)*phi1+(diff(y,lambda)-tau)*phi2+mu*diff(y,phi1)+lambda*diff(y,phi2)
join(flatten entries hxEqns, flatten entries hyEqns)

polyToProj = (p1, p2) -> (
    S := (coefficientRing ring p1)[append((ring p1)_*,symbol u)];
    ideal(homogenize(sub(p1,S),u), homogenize(sub(p2,S),u), u^(first degree p1))
    )

polyToHB = (p1, p2) -> (
    syz gens polyToProj(p1,p2)
    )

polyToDual = (p1, p2) -> (
    
    )

randomPoly = (R, d) -> (
    k := coefficientRing R;
    sum apply(toList(0..d), i->sum apply(flatten entries basis(i,R), m -> random(k)*m))
    )

bezierCurve = (R,lsPts) -> (
    t := R_0;
    xCoords := lsPts / (i -> i_0);
    yCoords := lsPts / (i -> i_1);
    d  := length lsPts - 1;
    Bx := sum apply(toList(0..d), i->(1-t)^(d-i)*t^i*xCoords_i);
    By := sum apply(toList(0..d), i->(1-t)^(d-i)*t^i*yCoords_i);
    (Bx, By, 1)
    )



--Question: Show the determinant of the Vandermonde
--is nonzero when the a_i's are distinct

--the determinant is zero if and only if the columns
--are linearly dependent

--To say vectors v_1...v_n are linearly dependent
--e.g., v_1=c_2*v_2+...+c_n*v_n
--c_1*v_1+...+c_n*v_n = 0 and 
--*not all the c_i are zero*
--In our case the v_i are the column vectors

--Each row, say row j will be of the form
--p(a_j)=c_1*1+c_2*a_j+c_3*a_j^2+...+c_n*a_j^(n-1) = 0
--where p(x)=c_1*1+c_2*x+...+c_n*x^(n-1) = 0

--Claim: p(x) is the zero polynomial.
--Why?: If p(x) is the zero polynomial, then all
--the c_i are zero. Hence, the column vectors are
--*independent*, contrary to assumption of column
--vectors being *dependent*, i.e., det vanishing

--proof: The degree of p(x)=n-1. But each p(a_j) is
--is a root. There are n choices of j, so p(x) has 
--at least n roots! By <some cor> saying
--a degree d nonzero polynomial can have at most d 
--roots, p(x) is the zero polynomial. 

--Fun Idea: polynomial interpolation. Say I have
--points (x_1, y_1), ... ,(x_n,y_n) and I want to 
--find a one variable polynomial passing through them. With
--the condition that all the x_i are distinct.
--p(x)=(x-x_1)...(x-x_n)
--p(x_1)=0

--We are in the search of a polynomial of the form
--p(x)=c_0+c_1x+...+c_(n-1)x^(n-1)
--Claim: Through any n points there is a polynomial 
--of degree n-1 passing through those points
--p(x_1)=c_0+c_1x_1+...+c_(n-1)x_1^(n-1)
--...
--p(x_n)=c_0+c_1x_n+...+c_(n-1)x_n^(n-1)

--If we let V be the vandermonde matrix with a_i=x_i
--and we let c be the column vector [c_1 ... c_n]
--Vc=y where y is the column vector [y_1 ... y_n]

-- I = (x, y) is not a principal ideal.
-- Suppose x = fg for polynomials f and g. 
-- x | f or x | g. <--- justify with "x is a prime element"
-- Suppose x | f. Then h(x)x=f(x) for some poly h
-- x = h(x)xg(x), so 1=h(x)g(x). 
--But 0=deg(1)=deg(h(x))+deg(g(x)). Hence,
--deg(g(x))=deg(h(x))=0, so g(x) is a constant.  

-- Suppose x = fg for polynomials f and g.
-- Then deg(x) = 1, so deg(fg)=deg(f)+deg(g)=1
-- So deg(f) = 0 or deg(g) = 0. Hence, either f 
-- or g is a constant. 

--1.5.4. Let h = gcd(f,g). Claim: Af+Bg=h.

-- proof: By <prop 6> (h)=(f,g). In particular,
-- h is in (f,g), so h = Af+Bg for some A and B in k[x]

--When I use the word "effective" or "effective way"
--what I mean is there is an actual algorithm to compute it
--There is an effective way to find the A and B
--It's called the "Extended Euclidean Algorithm"

--Question: Is every ideal in k[x] principal?

--Question: How does (gcd(f,g)) relate to (f, g)?
-- (f,g)==(gcd(f,g)) 


--Example: Choose f(x)=x(x+1) g(x)=x. Then gcd(f,g)=x.
--h/f=x/x(x+1)=1/(x+1) isn't in k[x]. f/h = (x+1)
-- h/2f = 1/2(x+1) 

--1.5.5 If f and g are in k[x] then 
--(f,g) = (f-qg, g) for any q in k[x]

--proof: It suffices to show both inclusions
-- To show (f,g) <= (f-qg, g) suffices to show that
-- both f and g are in (f-qg, g). g is clearly in (f-qg,g)
-- Since g is in (f-qg,g), qg is in (f-qg,g) by
-- def of an ideal. So f=(f-qg)+qg is in (f-qg, g)


--Computes q and r such that f = q*g+r and 
--either r = 0 or deg(r) < deg(g). 
euclidDiv = (f, g) -> (
    R := ring g;
    q := 0;
    r := f;
    print(q, r);
    while ((r != 0) and (leadTerm(r) % leadTerm(g) == 0)) do (
	q = q + sub(leadTerm(r)/leadTerm(g), R);
	r = r - sub(leadTerm(r)/leadTerm(g), R)*g;
	print(q, r);
	);
    (q,r)
    )

--for the gcd, follow the algorithm in the book (I think
--in prop 6?) and then for gcd(f1,f2,...,fn), just iterate
--your gcd algorithm for two functions.

--What we are really doing by computing gcd(f1..fn) is
--finding the generator (up to constant) of the ideal 
--(f1...fn)

--so to check if a polynomial g is in an ideal I in k[x].
--I=(f1..fn)=(f). g is in I is the same thing as 
--g is a multiple of f. And so we just have to divide g by 
--f and get remainder zero. 


--This outputs (q, r) such that f=q*g+r 
-- and either r = 0 or deg(r)<deg(g)
euclidDiv = (f, g) -> (
    R := ring g;
    q := 0;
    r := f;
    i := 0;
    print(i,q,r);
    while ((r != 0) and (leadTerm(r) % leadTerm(g) == 0)) do (
	q = q + (leadTerm(r)//leadTerm(g));
	r = r - (leadTerm(r)//leadTerm(g))*g;
	i = i+1;
	print(i,q, r);
	);
    (q,r)
    )

--Question: If I write polynomial division as 

--       __________
--  A   |    B
--
-- Then is f A or is f B?
-- 

-- Call the initial assignment of q and r q_0 and r_0
-- Call the ith assigment after q_i and r_i


-- Question: In terms of the "diagram" for doing polynomial
-- long division what is q_1 and r_1?

-- Question: Why is the algorithm correct at each step?

-- Question: Why does the algorithm terminate?
-- Question: What were the key ingredients to the alg?






-- What should a monomial order have?

-- If I have a polynomial p = c_1*m_1 + ... + c_n*m_n
-- where m_i are monomials and c_i are scalars then 
-- m_1 > m_2 > ... > m_n under the monomial order

-- Question: What are properties such a > should have?
-- Goal: construct a monomial order

-- Properties: 1) A > B or A = B or A < B
-- 2) Transitivity: A > B and B > C then A > C 
-- 

-- Question: If two monomials m_i and m_j are distinct,
-- then should m_i = m_j under the ordering?

-- No. We want m_i = m_j if and only if m_i == m_j









--Question: How many monomial orders are there on k[x]?
-- Prove your answer using properties of monomial orderings


--Define X^alpha > X^beta if the last nonzero entry of 
--alpha-beta is negative. 

-- xyz >xyz^2  <--- corresponds to ---> (1,1,1) > (1,1,2)
-- (1,1,1)-(1,1,2) = (0,0,-1)

--pure reverse lexicographic ordering

-- x^2y > xy <--> (2,1) > (1,1) FALSE
-- (2,1) - (1,1) = (1,0)

-- xy > x^2y

-- xy^2 < x^2y

-- You are minimizing bad x_1 > ... > x_n

-- xy > xy^2 > xy^3 > xy^4 > ...

-- 1 > x > x^2 > ...

-- Question: Think about why "graded reverse lexicographic ordering"
-- DOES give you a monomial order




--Gradings.
--ZZ-grading. 

--Polynomial Ring S = QQ[x_1..x_5]

--To get a ZZ-graded structure, I should figure out
--which graded piece every element is in. With that in 
--mind, how can I assign an integer to every polynomial?

--Definition: Let p be a polynomial. The support of p, 
--denoted supp(p), is the set of all monomials with nonzero
--coefficient in p.

--DEFINITION: A polynomial is homogeneous if all of the 
--monomials in its *support* are the same degree. 

--DEGREE: Let S_i be the set of all *homogeneous* polynomials
--of degree i. 

--Example: QQ[x,y]
--Question: Is x^2 in S_2 as defined above? Yes
--Question: Is -x^2+x in S_2 as defined above? No

--Question: Is S_i an additive group?
-- Yes. In fact, S_i is a QQ-vector space with basis 
-- given by the monomials of degree i.

--Definition: Define the *hilbert function* of S = kk[x_1..x_n]
--to be the function f : Z -> Z given by 
-- f(i) = dim_kk S_i

--f(i) = # monomials of degree i in S

--  A degree d monomial could be specified by d stars and 
--  n-1 bars.
--  i.e, a degree 5 monomial by **||** in kk[x,y,z]
--  So there are d + (n-1) total slots and then we just
--  need to choose either the d stars or the n-1 bars.

--Hence, f(i) = binomial(d+(n-1), d) = binomial(d+ (n-1), n-1)

--In general, let R = kk[x_1..x_n] / I
-- Then the hilbert function of R, denoted HF_R(i) is 
--the kk-vector space dimension of [R]_i, where 

--[R]_i is just the image of [S]_i under the projection map
-- S = k[x_1..x_n] --> k[x_1..x_n]/I = R

--[S]_2 = <x^2, x*y, y^2>_QQ
--[R]_2 = <x^2, y^2>_QQ

--Sloppiness: Why does the projection map preserve grading?
--Answer: It doesn't always, but if I is a homogeneous ideal, then 
--it does.

--Question: Suppose I have two rings [R]_i and [S]_i
--If HF(R,i) = HF(S,i) for all i then do they equal?

--Well, as sets they might be distinct, for example
--maybe the base fields of R and S are different but the
--dimensions are the same.

--Question: " Isomorphic?

--A finitely generated k-algebra is just a quotient
--of a polynomial ring over k.

--Setup: R = k[x_1..x_n]/I and I don't know what
-- I is. But I know certain relations are satisfied
-- in R. That means, I know the ideal J generated by
-- these relations is contained in I.
-- I also know the Hilbert function of I!!!
-- Since J is contained in I, J == I if and only
-- if HF(J,i) == HF(I,i) for all i.
-- Moreover, by a theorem I proved you only 
-- need to check HF(J,3) == HF(I,3).

--R = k[t^2, t^3] = k[x,y]/ideal()










--By a ZZ-grading, what I mean is decomposing S into a direct
--sum of S_i such that 
--1) Each S_i is a QQ-vector space
--2) The intersect(S_i,S_j) = 0
--3) S = direct sum of the S_i
--4) The S_i are indexed by ZZ
--5) S_i*S_j subset S_(i+j)






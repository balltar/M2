--Parsing using: https://yiboyang.github.io/graphrel/
--Just input the string on the right into stringToIdeal below
getStringList = (string) -> (
    firstParse = separate("],",string);
    fp         = length firstParse-1;
    for i from 0 to fp list(
    	s = separate(":[",firstParse_i);
    	if i != fp then s_1 else (separate("]",s_1))_0
    )
)

stringToIdeal = (string) -> (
    --parse string
    ls = getStringList(string);
    V  = QQ[vars (0..length ls - 1)];
    --print vars V;
    
    
    --make into ideal
    l = length ls-1;
    --print l;
    use V;
    print i;
    (V, ideal(flatten(for i from 0 to l list(
	lls  = flatten entries matrix ls_i;
	--print lls;
	ll   = length lls-1;
	--print ll;
	J = for j from 0 to ll list (V_i*V_(lls_j));
	--print J;
	--change next line to just J if don't want squares
	flatten {J,V_i*V_i}
	))))
)

isUnimodal = (testls) -> (
    l           = length testls - 1;
    i           = 0;
    goinguphill = true;
    while i < l do (
	if goinguphill then (
	    if (testls_i <= testls_(i+1)) then i=i+1
	    else goinguphill = false
	    )
	else if (testls_i >= testls_(i+1) and not goinguphill) then i=i+1
	else return false;
	);
    true
    )

isuhf = (nn,II) -> (
    M       = comodule II;
    testseq = for i from 0 to nn list hilbertFunction(i,M);
    (isUnimodal(testseq),testseq)
    )

processString = (nn, string) -> (
    (V, E) = stringToIdeal(string);
    isuhf(nn,E)
    )


testRingtest = ZZ/31[x,y,z,w];
testStringtest = "0:[1],1:[0,3,4],3:[1],4:[1]"; --star with center 1 and three branches 0,3,4
assert(stringToIdeal(testStringtest,testRingtest)==ideal(x*x,x*y,y*y,y*z,x*z,x*w))

--Generate
vertexRing = (n) -> (
    nn    = n-1;
    V = QQ[vars (0..nn)];
    V/ideal(for i from 0 to nn list V_i*V_i)
    )

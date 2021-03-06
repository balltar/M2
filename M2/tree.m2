restart
needsPackage "Nauty"
hilbertFunct = I -> (
    toList apply(0..15, i->hilbertFunction(i,comodule I))
    )

sameHBF = n -> (
    G := apply(generateGraphs(QQ[x_1..x_n],n-1), i->edgeIdeal i);
    P := apply(G, i->(hilbertFunct i,{i}));
    values hashTable(join, P)
    )

sameHBFTrees = n -> (
    G := apply(generateGraphs(QQ[x_1..x_n],n-1,OnlyConnected=>true), i->edgeIdeal i);
    P := apply(G, i->(hilbertFunct i,{i}));
    values hashTable(join, P)
    )

treeHilbFuncts = n -> (
    G := apply(generateGraphs(QQ[x_1..x_n],n-1,OnlyConnected=>true), i->edgeIdeal i);
    P := apply(G, i->(hilbertFunct i,{i}));
    fold(apply(keys hashTable(join, P),i-> net i), stack)
    )

needsPackage "Graphs"
--Given an edge ideal, output the graph
edgeIdealToGraph = I -> (
    ls := apply(I_*, i->{value (factor i)#0, value (factor i)#1});
    use ring I;
    --print (ring I === ring((ls_0)_0));
    graph((ring I)_*,ls)
    ) 

sameHBFGraphs = n -> (
    apply(sameHBF n, i->apply(i,j->edgeIdealToGraph j))
    )

--H#((keys H)#6)


needsPackage "Visualize"
visualizeList = (ls,name) -> (
    openPort "8080";
    p := "/home/tarball/" | name | "/";
    mkdir p;
    apply(ls, i -> visualize(i,VisPath => p, Warning => false));
    s = "rm -rf " | p | "*/";
    run s;
    closePort();
)

--generate just edge ideals of all trees on n vertices
generateTreesE = n -> (
    local i;
    S := QQ[symbol x_1..symbol x_n];
    apply(generateGraphs(S, n-1, OnlyConnected=>true), i->edgeIdeal(i))
    )

generateTreesM = n -> (
    local i;
    S := QQ[symbol x_1..symbol x_n];
    apply(generateGraphs(S,n-1,OnlyConnected=>true), i->intersect(edgeIdeal(i),(monomialIdeal(S_*))^[2]))
    )

--hvec
pvec = I -> (
    F = res I;
    for i from 0 to length F - 1 list rank(F#i)
    )

--hpoly
hpoly = I -> (
    numerator hilbertSeries I
    )

--Checks if a sequence S is peak unimodal
isUnimodal = S -> (
    l := #S-1;
    m := max(select(1..l, i -> S#(i-1)<S#i));
    M := min(select(1..l, i->S#(i-1)>S#i));
    f := M-m;
    if f<0 then all(f, j -> S#(m+1)==S#(m+f)) else true
    )

--Outputs a pair, one being the peak, other being length
peakIndex = S -> (
    l := #S-1;
    M := min(select(1..l, i-> S#(i-1)>=S#i));
    print S;
    (M,l+1,(M)/(l+1))
    )

checkMidpoints = n -> (
    els := generateTrees(n);
    apply(els, i-> peakIndex(apply(0..n, j->hilbertFunction(j, comodule i))))
    )

--Checks if a sequence S is log concave
isLogConcave = S -> (
    all(1..#S-2,i -> i*i>= (i-1)*(i+1))
    )

--Turns an integer sequence to a polynomial
seqToPoly = seq -> (
    R := ZZ[x];
    l := #seq-1;
    sum for i from 0 to l list seq_i*x^i   	
    )

postulationZero = (n,I) -> (
    P := hilbertPolynomial(I,Projective=>false);
    i := (ring(ideal P))_0;
    lsP := apply(n..10,j->sub(P,{i=>j}));
    lsH := apply(n..10,j->hilbertFunction(j,comodule I));
    print lsP;
    print lsH;
    lsP == lsH 
    )

--Generate algebra of perfect tree
perfTreeAlg = l -> (
    --Generate all 0-1 sequences of a certain length.
    generateIndices = l -> (
    prev  := {sequence (0)};
    union := prev;
    for i from 0 to l-1 do (
	next := flatten apply(prev, v->{join(v,sequence(0)),join(v,sequence(1))});
	prev = next;
	union = join(union, next);
	);
    union
    );

    indexes = generateIndices l;
    R = QQ[apply(indexes, i->x_i)];
    
    --Generate all relations in perfect binary tree
    generateRels = (l,R) -> (
    prev := {sequence(0)};
    rels := {};
    for i from 0 to l-1 do (
	use R;
	next := flatten apply(apply(prev,k->x_k), v->{v*x_(join((baseName v)#1,sequence 0)),v*x_(join((baseName v)#1,sequence 1))});
	prev = flatten apply(prev, v->{join(v,sequence(0)),join(v,sequence(1))});
	rels = join(rels, next);
	);
    monomialIdeal rels
    );

    I = generateRels(l,R);
    m = monomialIdeal R_*;
    I+(m^[2])
    )

--Generate algebra of perfect tree without the squares
perfTreeAlgEdge = l -> (
    indexes = generateIndices l;
    R = QQ[apply(indexes, i->x_i)];
    generateRels(l,R)
    )



needsPackage "Nauty"
--Generates all trees on n vertices with square of variables
generateTrees = n -> (
    local i;
    S = QQ[symbol x_1..symbol x_n];
    apply(generateGraphs(S, n-1, OnlyConnected=>true), i->edgeIdeal(i)+((ideal S_*)^[2]))
    )

--Gives map from [M]_s -> [M]_(s+1)
--So to compute [R/I]_s -> [R/I]_(s+1) use
--mapGradedPiece(s,comodule I)  
--Also, notice that the shifts are technically off, but since
--this is used for rank computation doesn't matter. If you
--wanted to fix this use something like:
--map(R^{(rank target mat):s+1},R^{(rank source mat):s}, mat),
--where mat is matrix constructed below and R=ring M
mapGradedPiece = (s,M) -> (
    ((sum (ring(M))_*)*basis(s,M))//basis(s+1,M)
    )

--Check WLP
WLPCheck = I -> (
    --Given a linear form l, a Z-graded module M, and a degree s, constructs the matrix
    --representing the linear map M_s -> M_(s+1) given by multiplication by l
    gradedPiece = 1;
    while(not (basis(gradedPiece,comodule I) == 0)) do (
	print "matrix construction:";
	M := time mapGradedPiece(gradedPiece,comodule I);
	print "rank:";
	r := time rank M;
	s := numgens source M;
	t := numgens target M;
	if (not (r==min(s, t))) then 
	    (return (false,gradedPiece));--,I));--(r,(t,s),false));--print(gradedPiece,r,(t,s),M,I);return false);
	gradedPiece = gradedPiece + 1;
	);
    print "\n";
    (true-1)--,I)
    )

delete(,apply(apply(generateTrees 8, j -> WLPCheck j), i-> if not i_0 then (if not hasZeroCols mapGradedPiece(i_1,comodule i_2) then (i_1,i_2))))
apply(delete(,apply(apply(generateTrees 8, j -> WLPCheck j), i-> if not i_0 then (if not hasZeroCols mapGradedPiece(i_1,comodule i_2) then (i_1,i_2))))
,i->mapGradedPiece(i_0,comodule i_1))
apply(apply(delete(,apply(apply(generateTrees 8, j -> WLPCheck j), i-> if not i_0 then (if not hasZeroCols mapGradedPiece(i_1,comodule i_2) then (i_1,i_2)))), i-> mapGradedPiece(i_0,comodule i_1)),j -> rank source j > rank target j)

hasZeroCols = mat -> (
    not select(apply(0..rank source mat-1, i->mat_i), j->zero sum flatten entries j) == ()
    );

pathIdeal = n -> (
    R := QQ[x_1..x_n];
    (ideal(R_*))^[2]+ideal toList apply(1..n-1, i->x_i*x_(i+1))
    )

starIdeal = n -> (
    R := QQ[x_1..x_n];
    (ideal(R_*))^[2]+ideal toList apply(2..n, i-> x_i*x_1)
    )
--roots hvec(J)
--binomRow = n -> for i from 0 to n list binomial(n,i)
    
hlength = I -> (
    count := 0;
    while true do (
	if hilbertFunction(count,comodule I) == 0 then break
	else count = count+1;
	);
    return count;
    )

constructIndepPoset = I -> (
    l := hlength I;
    for i in 0..l-1 do (
	
	)
    )
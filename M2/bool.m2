--Input: A positve integer n
--Output: All binary sequences of length n
binarySeqs = n -> (
    binSeqRec(n, {})
    )

binSeqRec = (n, ls) -> (
    if zero n 
    then {ls} 
    else join(binSeqRec(n-1,append(ls, 0)), 
	binSeqRec(n-1,append(ls,1)))
    )

--Input: A positive integer n
--Output: The free boolean Z/2Z-algebra wit n gens
boolAlg = n -> (
    R := ZZ/2[symbol a_1..symbol a_n];
    local i;
    R / ideal((1..n) / (i -> a_i^2-a_i))
    )

vee = (x,y) -> (
    x+x*y+y
    )

neg = x -> 1+x

princ = I -> fold(I_*, vee)

edgeIdealToAnnSeq = I -> (
    gls := I_*;
    B := ring I;
    apply(gls, g->(g' := sub(g,B); B = B/ideal g'; princ ann g'))
    )
    

boolAlgElements = B -> (
    n := length B_*;
    local i,j;
    mons := join apply(0..n, i->flatten entries basis(i, B));
    apply(binarySeqs(length mons), i->sum(apply(toList(0..length i-1), j->i_j*mons_j)))
    )

orMultTable = B -> (
    elms = boolAlgElements B;
    local i,j;
    netList table(elms, elms, (i,j)->i+i*j+j)
    )

stateDescriptions = B -> apply(binarySeqs(length B_*), i->product apply(toList(0..length i-1), j->i_j*B_j+(1+i_j)*neg(B_j)))

secondLayer = sds -> apply(subsets(sds, 2), i->vee(i_0,i_1))

nextLayer = (gls,Bgens) -> (
    delete(, unique apply(subsets(gls, 2), i -> (
	    if member(i_0*(i_1+1), Bgens)
	    then vee(i_0,i_1)
	    )
	))
    )

stateDescriptionsRels = B -> (
    gls := delete(0_B, apply(binarySeqs(length B_*), i->product apply(toList(0..length i-1), j->i_j*B_j+(1+i_j)*neg(B_j))));
    rels := apply(gls, i-> {0_B, i});
    (gls, rels, gls)
    )

secondLayerRels = (sds, rels, sds) -> (
    newRels := rels;
    newGls := apply(subsets(sds, 2), i-> (
	    v := vee(i_0,i_1); 
	    newRels = append(newRels, {i_0,v});
	    newRels = append(newRels, {i_1,v});
	    v
	    )
	);
    (newGls, newRels, sds)
    )

nextLayerRels = (gls, rels, sds) -> (
    newRels := rels;
    newGls := delete(, unique apply(subsets(gls, 2), i -> (
	    if member(i_0*(i_1+1), sds)
	    then (
		v := vee(i_0,i_1);
		newRels = append(newRels, {i_0, v});
		newRels = append(newRels, {i_1,v});
		v
		)
	    )
	));
    (newGls, newRels, sds)
    )

getHasse = n -> (
    B = boolAlg n;
    sdsR := stateDescriptionsRels B;
    << "Layer 1 computed. Size is " << length sdsR_0 << endl;
    prevLayer := secondLayerRels sdsR;
    << "Layer 2 computed. Size is " << length prevLayer_0 << endl;
    i := 3;
    numGens := length sdsR_0;
    while length prevLayer_0 != numGens do (
	prevLayer = nextLayerRels prevLayer;
	<< "Layer " << i << " computed. Size is " << length prevLayer_0 << endl;
	i = i+1;
	);
    join(prevLayer_1, apply(prevLayer_0, i->{i,1_B})) 
    )



needsPackage "Graphs";
needsPackage "Visualize";

openPort "8080";
    

--Compute DFS of graph

genLattice = B -> (
    elms = boolAlgElements B;
    --print elms;
    rels = toList set flatten apply(elms, i -> apply(elms, j-> if i+i*j+j == j then {i,j} else {i,i}));
    G = digraph(rels);
    visualize transitiveReduction G
    )

closePort();

transitiveReduction = G -> (
    remainingVertices := vertexSet G;
    newGraph := G;
    --print newGraph;
    while not zero length remainingVertices do (
	v := first remainingVertices;
	N := delete(v,toList children(newGraph,v));
	--print "%%%%%%%%%%%%%%%%%%%%%%%%%";
	--print v;
	--print N;
	while not zero length N do (
	    u := first N;
	    newGraph = digraph(toList(edges newGraph - set apply(delete(u,toList descendents(G,u)), u->{v,u})));
	    --print v;
	    --print u;
	    --print newGraph;
	    N = drop(N,1);
	    );
	remainingVertices = drop(remainingVertices, 1);
	);
    newGraph
    )

needsPackage "Nauty"
generateTreesE = n -> (
    local i;
    S := ZZ/2[symbol x_1..symbol x_n];
    I := ideal((1..n) / (i->x_i^2-x_i));
    apply(generateGraphs(S, n-1, OnlyConnected=>true), i->sub(edgeIdeal(i),S/I))
    )

generateGraphsE = n -> (
    local i;
    S := ZZ/2[symbol x_1..symbol x_n];
    I := ideal((1..n) / (i->x_i^2-x_i));
    apply(generateGraphs(S, OnlyConnected=>true), i->sub(edgeIdeal(i),S/I))
    )


generateTreeGens = n -> apply(generateTreesE n, i->ideal fold(i_*,vee))
generateGraphsGens = n -> apply(generateGraphsE n, i->ideal fold(i_*,vee))
generateNonTreesGens = n -> (
    Gs := generateGraphsGens n;
    R := ring first Gs;
    toList(Gs - set (apply(generateTreeGens n, i->sub(i,R))))
    )

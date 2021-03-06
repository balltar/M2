--We construct a monomial resolution of a 
--Gorenstein-linear ideal. The way we do this is by
--constructing a mapping cone of the resolution of the 
--homogeneous maximal ideal and the resolution of the 
--canonical module (shifted by one degree) with the map
--between them given by the one induced by the macaulay inverse system

--Generate all strictly increasing p-tuples of integers
--between 1 and d, where d>=1 and p<=d.
generateAllIncreasing = {ls => (), start => 1} >> opt -> (p,d) -> (
    if zero p then return opt.ls else
    flatten for i from opt.start to d list generateAllIncreasing(start=>i+1,p-1,d,ls=>append(opt.ls,i))
)

--Generate all nondecreasing p-tuples of integers
--between 1 and d, where d>=1 and p<=d.
generateAllNonDecreasing = {ls => (), start => 1} >> opt -> (p,d) -> (
    if zero p then return opt.ls else
    flatten for i from opt.start to d list generateAllIncreasing(start=>i,p-1,d,ls=>append(opt.ls,i))
)

--given a (p,q) returns basis for L(p,q) in the form of ordered
--pairs, with the first entry a sequence of a's and the second
--entry a sequence of b's.
--NOTE: We allow for an empty a sequence. This case is denoted by (0).
genBasisL = (p,q,d) -> (
    lsA := generateAllIncreasing(p+1,d);
    lA  := length lsA;
    if zero lA then (lsA = {toSequence {0}}; lA =1;);
    flatten for i from 0 to lA-1 list (
	lsB := generateAllNonDecreasing(q-1,d,start=>lsA_i_0);
	lB  := length lsB;
	for j from 0 to lB-1 list (lsA_i,lsB_j)
	)
    )

--Given a (p,q) returns basis for K(p,q) in the form of ordered
--pairs, with the first entry a sequence of a's and the second 
--entry a sequence of b's.
--NOTE: We allow for an empty a sequence. This case is denoted by (0).
genBasisK = (p,q,d) -> (
    lsB := generateAllNonDecreasing(q+1,d);
    lB  := length lsB;
    print lsB;
    flatten for i from 0 to lB-1 list (
	lsA := generateAllIncreasing(d-p-1,d,start=>lsB_i_0+1);
	lA  := length lsA;
	if zero lA then (lsA = {toSequence {0}}; lA =1;);
	for j from 0 to lA-1 list (lsA_j,lsB_i)
	)
    )

mu = I -> (
    rank source mingens I
    )

randomGenerators = (embdim, starting, ending) -> (
    range = random new List from (starting..ending);
    ls := {range_0};
    for i from 2 to embdim do (
	for j in range do (
	    	if (not member(j,ls)) and gcd(append(ls,j)) == 1 then (ls = append(ls,j); break;);
	    );
	if (length ls)<i then error "No generating set of that embedding dimension exists in specified range";
	);
    print ls;
    ls
    )

defIdeal = ls -> (
    local i;
    l := length ls;
    R := QQ[vars (0..l-1), Degrees=>ls]; -- Weights => {1,1,1,1,1}];
    print degree R_0;
    x := R_0;
    S := QQ[t];
    sub(substitute(ideal flatten entries mingens ker map(S,R,matrix {for i from 0 to l-1 list t^(ls_i)}), (R/x)),QQ[delete(x,R_*),Degrees => delete(ls_0,ls)])
    --R' := R/(ideal(x));
    --I' := substitute(I,R');
    )

ginMod = {n => 0} >> opt -> (I) -> (
    S := ring I;
    R := coefficientRing S;
    vls := remove(S_*, opt.n);
    print vls;
    gin substitute(substitute(I,(S/S_(opt.n))),R[vls])
    ) 

initialDegree = I -> (
    i := 1;
    m := ideal((ring I)_*);
    while isSubset(I,m^i) and isSubset(I,m^(i+1)) do
    i = i + 1;
    i
    )

--Examples of symmetric gen 5 with |mingens|=13: 
--{19,23,29,31,37},{19,27,28,31,32},{23,28,32,45,54},{19,27,29,32,33}

isSymmetric = I -> (
    F := res I;
    print F;
    if rank(F#(length((ring I)_*)-1)) == 1 then true else false
    )

delete(0,(for i from 0 to 100 list (
    rr = randomGenerators(5,6,70);
    II = defIdeal rr;
    if isSymmetric II then (mu II,rr) else (0)
    )))

    
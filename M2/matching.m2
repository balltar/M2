needsPackage "SimplicialComplexes"
needsPackage "Depth"
---------Builds Perfect Matching
perfectMatching = n -> (
    R := QQ[m_1..m_n];
    simplicialComplex({product R_*})
    )

--given simplicial complex and list of adjacent faces given a symbol name sym
addEdge = (S,sym,Als) -> (
    Fls := flatten entries facets(S);
    R  := (ring Fls_0);
    m  := dim R;
    R' := coefficientRing(R)[R_*,sym];
    V := R'_(length R'_*-1);
    
    simplicialComplex(delete(0,flatten(for i from 0 to m-1 list (
		Fls' := for F in flatten entries faces(i,S) list sub(F,R');
    		Als' := for A in Als list sub(A,R');
    		join(Fls',{V},apply(Fls', F->if any(Als', A -> (F%A)==0) then 0 else F*V))
		))))
    )

-----------Examples
S = addEdge(addEdge(perfectMatching(3), f_1, {m_1}),e_1,{m_1,m_2})
S' = addEdge(addEdge(S,f_2,{m_2,e_1}),e_2,{m_2,m_3})
S'' = addEdge(addEdge(S',f_3,{m_2,e_2,m_3}),e_3,{f_2,m_2,m_3,e_1})

T = addEdge(addEdge(S'', m_4, {}),e_4,{m_3,e_3})
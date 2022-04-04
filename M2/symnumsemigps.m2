homogenizeDeg = (ls,u) -> (
    R := coefficientRing(ring ls_0)[append((ring ls_0)_*,u)];
    u2 := last R_*;
    hls := ls / (i->homogenize(sub(i,R),u2));
    m := max apply(hls, i->first degree i);
    apply(hls, i->i*u2^(m-first degree i))
    )

--semigroup gens to map
sgtm = (ls, vls) -> (
    try (
	t := first vls;
	u := last vls;
	local R;
	try (
	    R = coefficientRing(ring t)[t];
	    ideal homogenizeDeg(apply(ls, i->(sub(t,R))^i), u)
	    )
	else (
	    R = QQ[t];
	    ideal homogenizeDeg(apply(ls, i->(sub(t,R))^i), u)
	    )
	)
    else (
	t := vls;
	try (
	    R = coefficientRing(ring t)[t];
	    ideal apply(ls, i->(sub(t,R))^i)
	    )
	else (
	    R = QQ[t];
	    ideal apply(ls, i->(sub(t,R))^i)
	    )
	
	)
    )

specialFiberEqns = I -> (
    n := length I_*;
    R := QQ[T_0..T_(n-1), Degrees => apply(toList(0..n-1), i->first degree I_i)];
    ker map(ring I, R, gens I)
    )

defIdeal = ls -> specialFiberEqns sgtm(ls, t)

bsgtm = (ls, vls) -> (
    betti res specialFiberEqns sgtm(ls, vls)
    )

rsftm = (ls, vls) -> (
    res specialFiberEqns sgtm(ls,vls)
    )
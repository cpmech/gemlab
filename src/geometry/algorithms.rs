/*
// Computes the distance from point to line passing through a -> b
pub fn dist_point_line(p, a, b *Point, tol float64, verbose bool) float64 {
    ns := NewSegment(a, b)
    vs := NewSegment(p, a)
    nn := ns.Len()
    if nn < tol { // point-point distance
        if verbose {
            io.Pfred("basicgeom.go: DistPointLine: __WARNING__ point-point distance too small:\n p=%v a=%v b=%v\n", p, a, b)
        }
        return vs.Len()
    }
    n := ns.Vector(1.0 / nn)
    v := vs.Vector(1.0)
    s := VecDot(v, n)
    l := VecNewAdd(1, v, -s, n) // l := v - dot(v,n) * n
    return VecNorm(l)
}
*/

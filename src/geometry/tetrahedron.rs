/// Calculates the tetrahedron (four) coordinates
///
/// # Output
///
/// * `zeta` -- the tetrahedron coordinates (len = 4)
///
/// # Input
///
/// * `xa,xb,xc,xd` -- corners (each has len = 3)
/// * `xp` -- the point to calculate the coordinates (len = 3)
///
/// # Panics
///
/// This function will panic if the array sizes are incorrect
pub fn tetrahedron_coords(zeta: &mut [f64], xa: &[f64], xb: &[f64], xc: &[f64], xd: &[f64], xp: &[f64]) {
    let (x1, y1, z1) = (xa[0], xa[1], xa[2]);
    let (x2, y2, z2) = (xb[0], xb[1], xb[2]);
    let (x3, y3, z3) = (xc[0], xc[1], xc[2]);
    let (x4, y4, z4) = (xd[0], xd[1], xd[2]);
    let (x, y, z) = (xp[0], xp[1], xp[2]);
    let v6 = (x2 - x1) * ((y2 - y3) * (z3 - z4) - (y3 - y4) * (z2 - z3))
        + (x3 - x2) * ((y3 - y4) * (z1 - z2) - (y1 - y2) * (z3 - z4))
        + (x4 - x3) * ((y1 - y2) * (z2 - z3) - (y2 - y3) * (z1 - z2));
    zeta[0] = ((x2 * (y3 * z4 - y4 * z3) + x3 * (y4 * z2 - y2 * z4) + x4 * (y2 * z3 - y3 * z2))
        + x * ((y4 - y2) * (z3 - z2) - (y3 - y2) * (z4 - z2))
        + y * ((x3 - x2) * (z4 - z2) - (x4 - x2) * (z3 - z2))
        + z * ((x4 - x2) * (y3 - y2) - (x3 - x2) * (y4 - y2)))
        / v6;
    zeta[1] = ((x1 * (y4 * z3 - y3 * z4) + x3 * (y1 * z4 - y4 * z1) + x4 * (y3 * z1 - y1 * z3))
        + x * ((y3 - y1) * (z4 - z3) - (y3 - y4) * (z1 - z3))
        + y * ((x4 - x3) * (z3 - z1) - (x1 - x3) * (z3 - z4))
        + z * ((x3 - x1) * (y4 - y3) - (x3 - x4) * (y1 - y3)))
        / v6;
    zeta[2] = ((x1 * (y2 * z4 - y4 * z2) + x2 * (y4 * z1 - y1 * z4) + x4 * (y1 * z2 - y2 * z1))
        + x * ((y2 - y4) * (z1 - z4) - (y1 - y4) * (z2 - z4))
        + y * ((x1 - x4) * (z2 - z4) - (x2 - x4) * (z1 - z4))
        + z * ((x2 - x4) * (y1 - y4) - (x1 - x4) * (y2 - y4)))
        / v6;
    zeta[3] = ((x1 * (y3 * z2 - y2 * z3) + x2 * (y1 * z3 - y3 * z1) + x3 * (y2 * z1 - y1 * z2))
        + x * ((y1 - y3) * (z2 - z1) - (y1 - y2) * (z3 - z1))
        + y * ((x2 - x1) * (z1 - z3) - (x3 - x1) * (z1 - z2))
        + z * ((x1 - x3) * (y2 - y1) - (x1 - x2) * (y3 - y1)))
        / v6;
}

/// Indicates if a point is inside a tetrahedron by looking at its tetrahedron coordinates (zeta)
///
/// Note: the point is inside (or on a face) if all zeta are positive (or zero)
///
/// # Input
///
/// * `zeta` -- the tetrahedron coordinates (len = 4)
///
/// # Panics
///
/// This function will panic if the array sizes are incorrect
#[inline]
pub fn in_tetrahedron(zeta: &[f64]) -> bool {
    if zeta[0] < 0.0 || zeta[1] < 0.0 || zeta[2] < 0.0 || zeta[3] < 0.0 {
        false
    } else {
        true
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{in_tetrahedron, tetrahedron_coords};
    use crate::util::ONE_BY_3;
    use russell_chk::vec_approx_eq;

    #[test]
    fn tetrahedron_coords_works() {
        let xa = &[0.0, 0.0, 0.0];
        let xb = &[1.0, 0.0, 0.0];
        let xc = &[0.0, 1.0, 0.0];
        let xd = &[0.0, 0.0, 1.0];
        let mut zeta = vec![0.0; 4];
        tetrahedron_coords(&mut zeta, xa, xb, xc, xd, &[0.0, 0.0, 0.0]);
        assert_eq!(zeta, &[1.0, 0.0, 0.0, 0.0]);
        tetrahedron_coords(&mut zeta, xa, xb, xc, xd, &[1.0, 0.0, 0.0]);
        assert_eq!(zeta, &[0.0, 1.0, 0.0, 0.0]);
        tetrahedron_coords(&mut zeta, xa, xb, xc, xd, &[0.0, 1.0, 0.0]);
        assert_eq!(zeta, &[0.0, 0.0, 1.0, 0.0]);
        tetrahedron_coords(&mut zeta, xa, xb, xc, xd, &[0.0, 0.0, 1.0]);
        assert_eq!(zeta, &[0.0, 0.0, 0.0, 1.0]);

        tetrahedron_coords(&mut zeta, xa, xb, xc, xd, &[0.5, 0.5, 0.0]);
        assert_eq!(zeta, &[0.0, 0.5, 0.5, 0.0]);
        tetrahedron_coords(&mut zeta, xa, xb, xc, xd, &[0.5, 0.5, 0.5]);
        assert_eq!(zeta, &[-0.5, 0.5, 0.5, 0.5]);
        tetrahedron_coords(&mut zeta, xa, xb, xc, xd, &[1.0, 1.0, 1.0]);
        assert_eq!(zeta, &[-2.0, 1.0, 1.0, 1.0]);
        tetrahedron_coords(&mut zeta, xa, xb, xc, xd, &[ONE_BY_3, ONE_BY_3, ONE_BY_3]);
        vec_approx_eq(&zeta, &[0.0, ONE_BY_3, ONE_BY_3, ONE_BY_3], 1e-15);
    }

    #[test]
    fn in_tetrahedron_works() {
        assert_eq!(in_tetrahedron(&[-1.0, 0.0, 0.0, 0.0]), false);
        assert_eq!(in_tetrahedron(&[0.0, -1.0, 0.0, 0.0]), false);
        assert_eq!(in_tetrahedron(&[0.0, 0.0, -1.0, 0.0]), false);
        assert_eq!(in_tetrahedron(&[0.0, 0.0, 0.0, -1.0]), false);
        assert_eq!(in_tetrahedron(&[0.0, 0.0, 0.0, 0.0]), true);
    }
}

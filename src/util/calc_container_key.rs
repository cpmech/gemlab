/// Calculates the key of the container where the point should fall in
///
/// # Reference
///
/// * Durand, Farias, and Pedroso (2015) Computing intersections between
///   non-compatible curves and finite elements, Computational Mechanics;
///   DOI=10.1007/s00466-015-1181-y
#[inline]
pub fn calc_container_key(ndim: usize, side_length: f64, ndiv: &[usize], xmin: &[f64], x: &[f64]) -> usize {
    let mut ix = ((x[0] - xmin[0]) / side_length) as usize; // (Eq. 8)
    let mut iy = ((x[1] - xmin[1]) / side_length) as usize;
    if ix == ndiv[0] {
        ix -= 1; // point is on max edge => move to inner container
    }
    if iy == ndiv[1] {
        iy -= 1; // point is on max edge => move to inner container
    }
    if ndim == 2 {
        return ix + iy * ndiv[0];
    }
    let mut iz = ((x[2] - xmin[2]) / side_length) as usize;
    if iz == ndiv[2] {
        iz -= 1; // point is on max edge => move to inner container
    }
    return ix + iy * ndiv[0] + iz * ndiv[0] * ndiv[1];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::calc_container_key;

    #[test]
    fn calc_container_key_works() {
        // 2d
        let ndim = 2;
        let side_length = 0.30000000000000004;
        let ndiv = &[4, 8];
        let xmin = &[-0.30000000000000004, -0.30000000000000004];
        let x = &[0.1, 0.5];
        assert_eq!(calc_container_key(ndim, side_length, ndiv, xmin, x), 9);
        let x = &[0.7, 0.8];
        assert_eq!(calc_container_key(ndim, side_length, ndiv, xmin, x), 15);
        let x = &[-0.2, 1.8];
        assert_eq!(calc_container_key(ndim, side_length, ndiv, xmin, x), 24);
        let x = &[0.8, 1.8];
        assert_eq!(calc_container_key(ndim, side_length, ndiv, xmin, x), 27);
        let xmax = &[
            xmin[0] + (ndiv[0] as f64) * side_length,
            xmin[1] + (ndiv[1] as f64) * side_length,
        ];
        let x = &[xmax[0], xmax[1]];
        assert_eq!(
            calc_container_key(ndim, side_length, ndiv, xmin, x),
            (ndiv[0] * ndiv[1] - 1)
        );

        // 3d
        let ndim = 3;
        let side_length = 1.1;
        let ndiv = &[2, 2, 2];
        let xmin = &[-1.1, -1.1, -1.1];
        let x = &[-1.0, -1.0, -1.0];
        assert_eq!(calc_container_key(ndim, side_length, ndiv, xmin, x), 0);
        let x = &[1.0, 1.0, 1.0];
        assert_eq!(calc_container_key(ndim, side_length, ndiv, xmin, x), 7);
        let xmax = &[
            xmin[0] + (ndiv[0] as f64) * side_length,
            xmin[1] + (ndiv[1] as f64) * side_length,
            xmin[2] + (ndiv[2] as f64) * side_length,
        ];
        let x = &[xmax[0], xmax[1], xmax[2]];
        assert_eq!(
            calc_container_key(ndim, side_length, ndiv, xmin, x),
            (ndiv[0] * ndiv[1] * ndiv[2] - 1)
        );
    }
}

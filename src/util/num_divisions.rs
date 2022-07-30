use crate::StrError;

/// Computes the number of divisions for GridSearch aiming at a square/cubic containers
///
/// Note: **long** means the longest direction whereas
///       **other** corresponds to the `not-long` directions.
///
/// ```text
/// ndiv_other = truncate((delta_other/delta_long) * ndiv_long)
/// ndiv_other = max(ndiv_min, ndiv_other)
/// ```
///
/// # Input
///
/// * `ndiv_min` -- minimum number of divisions for the other directions (≥ 1)
/// * `ndiv_long` -- number of divisions for the longest direction (> ndiv_min)
/// * `xmin` -- min coordinates (len = ndim)
/// * `xmax` -- max coordinates (len = ndim); must be greater than min
///
/// # Output
///
/// * `ndiv` -- the number of divisions along each direction (ndim)
///
/// # Examples
///
/// ```
/// use gemlab::util::num_divisions;
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     let xmin = &[0.0, 0.0];
///     let xmax = &[1.0, 1.0];
///     let ndiv = num_divisions(2, 5, xmin, xmax)?;
///     assert_eq!(ndiv, &[5, 5]);
///
///     let xmax = &[1.0, 10.0];
///     let ndiv = num_divisions(2, 5, xmin, xmax)?;
///     assert_eq!(ndiv, &[2, 5]);
///     Ok(())
/// }
/// ```
pub fn num_divisions(ndiv_min: usize, ndiv_long: usize, xmin: &[f64], xmax: &[f64]) -> Result<Vec<usize>, StrError> {
    if ndiv_min < 1 {
        return Err("ndiv_min must be ≥ 1");
    }
    if ndiv_long <= ndiv_min {
        return Err("ndiv_long must be > ndiv_min");
    }
    if xmin.len() != xmax.len() {
        return Err("xmin.len() must be equal to xmax.len()");
    }
    let ndim = xmin.len();
    if ndim < 2 {
        return Err("ndim must be ≥ 2");
    }
    let delta: Vec<_> = xmin.iter().zip(xmax).map(|(a, b)| *b - *a).collect();
    let mut index_long = 0;
    let mut delta_long = delta[index_long];
    for i in 1..ndim {
        if delta[i] > delta_long {
            index_long = i;
            delta_long = delta[i];
        }
    }
    let mut ndiv = vec![0; ndim];
    for i in 0..ndim {
        if delta[i] <= 0.0 {
            return Err("xmax must be greater than xmin");
        }
        if i == index_long {
            ndiv[i] = ndiv_long;
        } else {
            ndiv[i] = ((delta[i] / delta_long) * (ndiv_long as f64)) as usize;
            ndiv[i] = usize::max(ndiv_min, ndiv[i]);
        }
    }
    Ok(ndiv)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::num_divisions;

    #[test]
    fn num_divisions_fails_on_errors() {
        assert_eq!(
            num_divisions(0, 20, &[0.0, 0.0], &[1.0, 1.0]).err(),
            Some("ndiv_min must be ≥ 1")
        );
        assert_eq!(
            num_divisions(2, 1, &[0.0, 0.0], &[1.0, 1.0]).err(),
            Some("ndiv_long must be > ndiv_min")
        );
        assert_eq!(
            num_divisions(2, 3, &[0.0, 0.0], &[1.0, 1.0, 1.0]).err(),
            Some("xmin.len() must be equal to xmax.len()")
        );
        assert_eq!(num_divisions(2, 3, &[0.0], &[1.0]).err(), Some("ndim must be ≥ 2"));
        assert_eq!(
            num_divisions(2, 20, &[2.0, 2.0], &[0.0, 0.0]).err(),
            Some("xmax must be greater than xmin")
        );
    }

    #[test]
    fn num_divisions_works() {
        assert_eq!(num_divisions(1, 10, &[0.0, 0.0], &[1.0, 1.0]).unwrap(), &[10, 10]);
        assert_eq!(num_divisions(2, 10, &[0.0, 0.0], &[1.0, 1.0]).unwrap(), &[10, 10]);
        assert_eq!(num_divisions(3, 20, &[0.0, 0.0], &[100.0, 100.0]).unwrap(), &[20, 20]);
        assert_eq!(num_divisions(3, 20, &[0.0, 0.0], &[1.0, 2.0]).unwrap(), &[10, 20]);
        assert_eq!(num_divisions(3, 20, &[0.0, 0.0], &[2.0, 1.0]).unwrap(), &[20, 10]);
        assert_eq!(num_divisions(3, 20, &[0.0, 0.0], &[1.0, 4.0]).unwrap(), &[5, 20]);
        assert_eq!(num_divisions(3, 20, &[0.0, 0.0], &[1.0, 5.0]).unwrap(), &[4, 20]);
        assert_eq!(num_divisions(2, 20, &[0.0, 0.0], &[1.0, 10.0]).unwrap(), &[2, 20]);
        assert_eq!(num_divisions(1, 20, &[0.0, 0.0], &[1.0, 100.0]).unwrap(), &[1, 20]);
        assert_eq!(num_divisions(2, 20, &[0.0, 0.0], &[1.0, 100.0]).unwrap(), &[2, 20]);
        assert_eq!(num_divisions(3, 20, &[0.0, 0.0], &[10000.0, 1.0]).unwrap(), &[20, 3]);
        assert_eq!(num_divisions(3, 20, &[0.0, 0.0], &[10.0, 3.0]).unwrap(), &[20, 6]);
        assert_eq!(num_divisions(3, 20, &[0.0, 0.0], &[10.0, 3.5]).unwrap(), &[20, 7]);
        assert_eq!(num_divisions(3, 20, &[-1.0, -2.0], &[1.0, 2.0]).unwrap(), &[10, 20]);
    }
}

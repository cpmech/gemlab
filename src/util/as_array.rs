/// Defines a trait to handle 1D arrays
///
/// # Example
///
/// ```
/// use gemlab::util::AsArray1D;
///
/// fn sum<'a, T, U>(array: &'a T) -> f64
/// where
///     T: AsArray1D<'a, U>,
///     U: 'a + Into<f64>,
/// {
///     let mut res = 0.0;
///     let m = array.size();
///     for i in 0..m {
///         res += array.at(i).into();
///     }
///     res
/// }
///
/// // heap-allocated 1D array (vector)
/// let x = vec![1.0, 2.0, 3.0];
/// assert_eq!(sum(&x), 6.0);
///
/// // heap-allocated 1D array (slice)
/// let y: &[f64] = &[10.0, 20.0, 30.0];
/// assert_eq!(sum(&y), 60.0);
///
/// // stack-allocated (fixed-size) 2D array
/// let z = [100.0, 200.0, 300.0];
/// assert_eq!(sum(&z), 600.0);
/// ```
pub trait AsArray1D<'a, U: 'a> {
    /// Returns the size of the array
    fn size(&self) -> usize;

    /// Returns the value at index i
    fn at(&self, i: usize) -> U;
}

/// Defines a heap-allocated 1D array (vector)
impl<'a, U: 'a> AsArray1D<'a, U> for Vec<U>
where
    U: 'a + Copy,
{
    fn size(&self) -> usize {
        self.len()
    }
    fn at(&self, i: usize) -> U {
        self[i]
    }
}

/// Defines a heap-allocated 1D array (slice)
impl<'a, U> AsArray1D<'a, U> for &'a [U]
where
    U: 'a + Copy,
{
    fn size(&self) -> usize {
        self.len()
    }
    fn at(&self, i: usize) -> U {
        self[i]
    }
}

/// Defines a stack-allocated (fixed-size) 1D array
impl<'a, U, const M: usize> AsArray1D<'a, U> for [U; M]
where
    U: 'a + Copy,
{
    fn size(&self) -> usize {
        self.len()
    }
    fn at(&self, i: usize) -> U {
        self[i]
    }
}

/// Defines a trait to handle 2D arrays
///
/// # Example
///
/// ```
/// use gemlab::util::AsArray2D;
///
/// fn sum<'a, T, U>(array: &'a T) -> f64
/// where
///     T: AsArray2D<'a, U>,
///     U: 'a + Into<f64>,
/// {
///     let mut res = 0.0;
///     let (m, n) = array.size();
///     for i in 0..m {
///         for j in 0..n {
///             res += array.at(i, j).into();
///         }
///     }
///     res
/// }
///
/// fn str_row<'a, T, U>(array: &'a T, i: usize) -> String
/// where
///     T: AsArray2D<'a, U>,
///     U: 'a + std::fmt::Debug,
/// {
///     format!("{:?}", array.row(i))
/// }
///
/// // heap-allocated 2D array (vector of vectors)
/// const IGNORED: f64 = 123.456;
/// let a = vec![
///     vec![1.0, 2.0],
///     vec![3.0, 4.0, IGNORED, IGNORED, IGNORED],
///     vec![5.0, 6.0],
/// ];
/// assert_eq!(sum(&a), 21.0);
/// assert_eq!(str_row(&a, 1), "[3.0, 4.0, 123.456, 123.456, 123.456]");
///
/// // heap-allocated 2D array (aka slice of slices)
/// let b: &[&[f64]] = &[
///     &[10.0, 20.0],
///     &[30.0, 40.0, IGNORED],
///     &[50.0, 60.0, IGNORED, IGNORED],
/// ];
/// assert_eq!(sum(&b), 210.0);
/// assert_eq!(str_row(&b, 0), "[10.0, 20.0]");
///
/// // stack-allocated (fixed-size) 2D array
/// let c = [
///     [100.0, 200.0],
///     [300.0, 400.0],
///     [500.0, 600.0],
/// ];
/// assert_eq!(sum(&c), 2100.0);
/// assert_eq!(str_row(&c, 0), "[100.0, 200.0]");
/// ```
pub trait AsArray2D<'a, U: 'a> {
    /// Returns the (m,n) size of the array
    ///
    /// # Panics
    ///
    /// This function panics if the array is empty.
    fn size(&self) -> (usize, usize);

    /// Returns the value at (i,j) indices
    ///
    /// # Panics
    ///
    /// This function panics if the indices are out of range.
    fn at(&self, i: usize, j: usize) -> U;

    /// Returns the i-th row
    ///
    /// # Panics
    ///
    /// This function panics if the index is out of range.
    fn row(&self, i: usize) -> &[U];
}

/// Defines a heap-allocated 2D array (vector of vectors)
///
/// # Notes
///
/// * The number of columns is defined by the first row
/// * The next rows must have at least the same number of columns as the first row
///
/// # Panics
///
/// The methods may panic if the array is empty.
impl<'a, U: 'a> AsArray2D<'a, U> for Vec<Vec<U>>
where
    U: 'a + Copy,
{
    fn size(&self) -> (usize, usize) {
        (self.len(), self[0].len())
    }
    fn at(&self, i: usize, j: usize) -> U {
        self[i][j]
    }
    fn row(&self, i: usize) -> &[U] {
        &self[i]
    }
}

/// Defines a heap-allocated 2D array (slice of slices)
///
/// # Notes
///
/// * The number of columns is defined by the first row
/// * The next rows must have at least the same number of columns as the first row
///
/// # Panics
///
/// The methods may panic if the array is empty.
impl<'a, U> AsArray2D<'a, U> for &'a [&'a [U]]
where
    U: 'a + Copy,
{
    fn size(&self) -> (usize, usize) {
        (self.len(), self[0].len())
    }
    fn at(&self, i: usize, j: usize) -> U {
        self[i][j]
    }
    fn row(&self, i: usize) -> &[U] {
        &self[i]
    }
}

/// Defines a stack-allocated (fixed-size) 2D array
///
/// # Panics
///
/// The methods may panic if the array is empty.
impl<'a, U, const M: usize, const N: usize> AsArray2D<'a, U> for [[U; N]; M]
where
    U: 'a + Copy,
{
    fn size(&self) -> (usize, usize) {
        (self.len(), self[0].len())
    }
    fn at(&self, i: usize, j: usize) -> U {
        self[i][j]
    }
    fn row(&self, i: usize) -> &[U] {
        &self[i]
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{AsArray1D, AsArray2D};
    use std::fmt::Write;

    fn array_1d_str<'a, T, U>(array: &'a T) -> String
    where
        T: AsArray1D<'a, U>,
        U: 'a + std::fmt::Display,
    {
        let mut buf = String::new();
        let m = array.size();
        for i in 0..m {
            write!(&mut buf, "{},", array.at(i)).unwrap();
        }
        write!(&mut buf, "\n").unwrap();
        buf
    }

    fn array_2d_str<'a, T, U>(array: &'a T) -> String
    where
        T: AsArray2D<'a, U>,
        U: 'a + std::fmt::Display,
    {
        let mut buf = String::new();
        let (m, n) = array.size();
        for i in 0..m {
            for j in 0..n {
                write!(&mut buf, "{},", array.at(i, j)).unwrap();
            }
            write!(&mut buf, "\n").unwrap();
        }
        buf
    }

    fn array_2d_str_row<'a, T, U>(array: &'a T, i: usize) -> String
    where
        T: AsArray2D<'a, U>,
        U: 'a + std::fmt::Debug,
    {
        format!("{:?}", array.row(i))
    }

    #[test]
    fn as_array_1d_works() {
        // heap-allocated 1D array (vector)
        let x_data = vec![1.0, 2.0, 3.0];
        assert_eq!(array_1d_str(&x_data), "1,2,3,\n");

        // heap-allocated 1D array (slice)
        let y_data: &[f64] = &[10.0, 20.0, 30.0];
        assert_eq!(array_1d_str(&y_data), "10,20,30,\n");

        // stack-allocated (fixed-size) 2D array
        let z_data = [100.0, 200.0, 300.0];
        assert_eq!(array_1d_str(&z_data), "100,200,300,\n");
    }

    #[test]
    fn as_array_2d_works() {
        // heap-allocated 2D array (vector of vectors)
        const IGNORED: f64 = 123.456;
        let a_data = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0, IGNORED, IGNORED, IGNORED],
            vec![5.0, 6.0],
        ];
        assert_eq!(
            array_2d_str(&a_data),
            "1,2,\n\
             3,4,\n\
             5,6,\n"
        );
        assert_eq!(array_2d_str_row(&a_data, 1), "[3.0, 4.0, 123.456, 123.456, 123.456]");

        // heap-allocated 2D array (aka slice of slices)
        let b_data: &[&[f64]] = &[&[10.0, 20.0], &[30.0, 40.0, IGNORED], &[50.0, 60.0, IGNORED, IGNORED]];
        assert_eq!(
            array_2d_str(&b_data),
            "10,20,\n\
             30,40,\n\
             50,60,\n"
        );
        assert_eq!(array_2d_str_row(&b_data, 0), "[10.0, 20.0]");

        // stack-allocated (fixed-size) 2D array
        let c_data = [[100.0, 200.0], [300.0, 400.0], [500.0, 600.0]];
        assert_eq!(
            array_2d_str(&c_data),
            "100,200,\n\
             300,400,\n\
             500,600,\n"
        );
        assert_eq!(array_2d_str_row(&c_data, 2), "[500.0, 600.0]");
    }
}

use super::*;
use russell_lab::Vector;

/// Defines the kind of shape
pub enum Kind {
    Lin2,
    Lin3,
    Lin5,
    Tri3,
    Tri6,
    Tri15,
    Qua4,
    Qua8,
    Qua9,
    Tet4,
    Tet10,
    Hex8,
    Hex20,
    Joint,
    Lin4,
    Tri10,
    Qua12,
    Qua16,
    Beam,
}

/// Defines the functionality of shape
pub trait Trait {
    /// Calculates the interpolation functions at natural coordinate ξ
    ///
    /// ```text
    /// interp[m](ξ) = Sᵐ(ξ)
    /// ```
    fn calc_interp(&mut self, coord: &Vector);

    /// Calculates the derivatives of interpolation fn at natural coordinate ξ
    ///
    /// ```text
    /// deriv[m][i](ξ) = ({dSᵐ(ξ)/dξ}_ξ)[i]
    /// ```
    fn calc_deriv(&mut self, coord: &Vector);

    /// Returns the previously computed interpolation fn for point m
    ///
    /// ```text
    /// interp[m](ξ) = Sᵐ(ξ)
    /// ```
    fn get_interp(&self, m: usize) -> f64;

    /// Returns the previously computed derivative of interpolation fn for point m
    ///
    /// ```text
    /// deriv[m][i](ξ) = ({dSᵐ(ξ)/dξ}_ξ)[i]
    /// ```
    fn get_deriv(&self, deriv: &mut Vector, m: usize);

    /// Returns the number of dimensions
    fn get_ndim(&self) -> usize;

    /// Returns the number of points
    fn get_npoint(&self) -> usize;

    /// Returns the number of edges
    fn get_nedge(&self) -> usize;

    /// Returns the number of faces
    fn get_nface(&self) -> usize;

    /// Returns all natural coordinates @ point m
    ///
    /// ```text
    /// coord = ξᵐ = vector{r, s, t} @ point m
    /// ```
    fn get_ksi(&self, ksi: &mut Vector, m: usize);
}

pub fn new(kind: Kind) -> Box<dyn Trait> {
    match kind {
        Kind::Qua4 => Box::new(Qua4::new()),
        Kind::Hex8 => Box::new(Hex8::new()),
        _ => panic!("Geo is not available yet"),
    }
}

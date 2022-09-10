use crate::shapes::Scratchpad;

use super::IntegPointData;

/// Holds common arguments for the numerical integration functions
pub struct CommonArgs<'a> {
    /// The temporary interpolation, Jacobian, gradient  variables computed at an integration point
    pub pad: &'a mut Scratchpad,

    /// Integration points (n_integ_point)
    pub ips: IntegPointData,

    /// Fills the output vector or matrix with zeros, otherwise accumulate values
    pub clear: bool,

    /// Specifies an axisymmetric model
    ///
    /// Performs the integration for 1 radian.
    pub axisymmetric: bool,

    /// Extra coefficient that can be used to define the thickness of plane-stress models
    pub alpha: f64,

    /// Stride marking the first row in the output vector/matrix where to add components
    pub ii0: usize,

    /// Stride marking the first column in the output matrix where to add components
    pub jj0: usize,
}

impl<'a> CommonArgs<'a> {
    /// Allocates new instance
    pub fn new(pad: &'a mut Scratchpad, ips: IntegPointData) -> Self {
        CommonArgs {
            pad,
            ips,
            clear: true,
            axisymmetric: false,
            alpha: 1.0,
            ii0: 0,
            jj0: 0,
        }
    }
}

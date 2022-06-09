use crate::shapes::{GeoClass, GeoKind};
use crate::StrError;

/// Defines an alias for integration points data (coordinates and weights)
///
/// Each integration point (IP) is defined by 3 reference
/// coordinates (r, s, t) and the weight (w).
///
/// The 1D and 2D geometries have `s` and/or `t` equal to zero as appropriate.
///
/// The data structure is a 2D array such that [[f64; 4]; NIP]
/// where `NIP` is the number of integration points in a particular set.
/// Therein, `4` corresponds to (r, s, t) and the weight (w).
pub type IntegPointData = &'static [[f64; 4]];

/// Returns a default set of integration points and weights
///
/// **Note:** The default set is chosen based on the number of nodes
/// and the "probable" use case in finite element analyses.
///
/// # Examples
///
/// ```
/// use gemlab::integ::default_integ_points;
/// use gemlab::shapes::GeoKind;
///
/// let ips = default_integ_points(GeoKind::Tet4);
/// assert_eq!(ips, [[0.25, 0.25, 0.25, 1.0/6.0]]);
/// ```
pub fn default_integ_points(kind: GeoKind) -> IntegPointData {
    match kind {
        // Lin
        GeoKind::Lin2 => &IP_LIN_LEGENDRE_2,
        GeoKind::Lin3 => &IP_LIN_LEGENDRE_3,
        GeoKind::Lin4 => &IP_LIN_LEGENDRE_4,
        GeoKind::Lin5 => &IP_LIN_LEGENDRE_5,
        // Tri
        GeoKind::Tri3 => &IP_TRI_INTERNAL_1,
        GeoKind::Tri6 => &IP_TRI_INTERNAL_4,
        GeoKind::Tri10 => &IP_TRI_INTERNAL_12,
        GeoKind::Tri15 => &IP_TRI_INTERNAL_16,
        // Qua
        GeoKind::Qua4 => &IP_QUA_LEGENDRE_4,
        GeoKind::Qua8 => &IP_QUA_LEGENDRE_9,
        GeoKind::Qua9 => &IP_QUA_LEGENDRE_9,
        GeoKind::Qua12 => &IP_QUA_LEGENDRE_16,
        GeoKind::Qua16 => &IP_QUA_LEGENDRE_16,
        GeoKind::Qua17 => &IP_QUA_LEGENDRE_16,
        // Tet
        GeoKind::Tet4 => &IP_TET_INTERNAL_1,
        GeoKind::Tet10 => &IP_TET_FELIPPA_14,
        GeoKind::Tet20 => &IP_TET_FELIPPA_24,
        // Hex
        GeoKind::Hex8 => &IP_HEX_LEGENDRE_8,
        GeoKind::Hex20 => &IP_HEX_LEGENDRE_27,
        GeoKind::Hex32 => &IP_HEX_LEGENDRE_27, // TODO: Implement higher order points for Hex32
    }
}

/// Selects integration points constants (coordinates and weights)
///
/// # Input
///
/// * `class` -- The geometry class
/// * `n_integ_point` -- Number of integration points desired (see Options below)
///
/// # Options
///
/// `n_integ_point` for **Lin** class:
///
/// * `1` -- Conventional Legendre integration points and weights
/// * `2` -- Conventional Legendre integration points and weights
/// * `3` -- Conventional Legendre integration points and weights
/// * `4` -- Conventional Legendre integration points and weights
/// * `5` -- Conventional Legendre integration points and weights
///
/// `n_integ_point` for **Tri** class:
///
/// * `1` -- Internal integration points and weights
/// * `3` -- Internal integration points and weights
/// * `1_003` -- Edge integration points and weights
/// * `4` -- Internal integration points and weights
/// * `12` -- Internal integration points and weights
/// * `16` -- Internal integration points and weights
///
/// `n_integ_point` for **Qua** class:
///
/// * `1` -- Conventional Legendre integration points and weights
/// * `4` -- Conventional Legendre integration points and weights
/// * `5` -- Wilson's integration points and weights. "Corner" version
/// * `1_005` -- 5 points. Wilson's integration points and weights. "Stable" version version with w0=0.004 and wa=0.999 to mimic 4-point rule
/// * `8` -- Wilson's integration points and weights.
/// * `9` -- Conventional Legendre integration points and weights
/// * `16` -- Conventional Legendre integration points and weights
///
/// `n_integ_point` for **Tet** class:
///
/// * `1` -- Internal integration points and weights, degree 1
/// * `4` -- Internal integration points and weights, degree 2 (based on Felippa's code; see bottom of source file)
/// * `5` -- Internal integration points and weights, degree 3 (with negative weight)
/// * `8` -- Felippa's ips and weights, degree 3 (based on Felippa's code; see bottom of source file)
/// * `14` -- Felippa's ips and weights, degree 4 (based on Felippa's code; see bottom of source file)
/// * `15` -- Felippa's ips and weights, degree 5 (based on Felippa's code; see bottom of source file)
/// * `24` -- Felippa's ips and weights, degree 6 (based on Felippa's code; see bottom of source file)
///
/// `n_integ_point` for **Hex** class:
///
/// * `6` -- Iron's integration points and weights
/// * `8` -- Conventional Legendre integration points and weights
/// * `9` -- Wilson's integration points and weights. "Corner" version
/// * `1_009` -- Wilson's integration points and weights. "Stable" version
/// * `14` -- Iron's integration points and weights
/// * `27` -- Conventional Legendre integration points and weights
///
/// # Examples
///
/// ```
/// use gemlab::integ::select_integ_points;
/// use gemlab::shapes::GeoClass;
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     let ips = select_integ_points(GeoClass::Tet, 1)?;
///     assert_eq!(ips, [[0.25, 0.25, 0.25, 1.0/6.0]]);
///     Ok(())
/// }
/// ```
pub fn select_integ_points(class: GeoClass, n_integ_point: usize) -> Result<IntegPointData, StrError> {
    let ips: IntegPointData = match class {
        // Lin
        GeoClass::Lin => match n_integ_point {
            1 => &IP_LIN_LEGENDRE_1,
            2 => &IP_LIN_LEGENDRE_2,
            3 => &IP_LIN_LEGENDRE_3,
            4 => &IP_LIN_LEGENDRE_4,
            5 => &IP_LIN_LEGENDRE_5,
            _ => return Err("desired number of integration points is not available for Lin class"),
        },
        // Tri
        GeoClass::Tri => match n_integ_point {
            1 => &IP_TRI_INTERNAL_1,
            3 => &IP_TRI_INTERNAL_3,
            1_003 => &IP_TRI_EDGE_3,
            4 => &IP_TRI_INTERNAL_4,
            12 => &IP_TRI_INTERNAL_12,
            16 => &IP_TRI_INTERNAL_16,
            _ => return Err("desired number of integration points is not available for Tri class"),
        },
        // Qua
        GeoClass::Qua => match n_integ_point {
            1 => &IP_QUA_LEGENDRE_1,
            4 => &IP_QUA_LEGENDRE_4,
            5 => &IP_QUA_WILSON_CORNER_5,
            1_005 => &IP_QUA_WILSON_STABLE_5,
            8 => &IP_QUA_WILSON_8,
            9 => &IP_QUA_LEGENDRE_9,
            16 => &IP_QUA_LEGENDRE_16,
            _ => return Err("desired number of integration points is not available for Qua class"),
        },
        // Tet
        GeoClass::Tet => match n_integ_point {
            1 => &IP_TET_INTERNAL_1,
            4 => &IP_TET_INTERNAL_4,
            5 => &IP_TET_INTERNAL_5,
            8 => &IP_TET_FELIPPA_8,
            14 => &IP_TET_FELIPPA_14,
            15 => &IP_TET_FELIPPA_15,
            24 => &IP_TET_FELIPPA_24,
            _ => return Err("desired number of integration points is not available for Tet class"),
        },
        // Hex
        GeoClass::Hex => match n_integ_point {
            6 => &IP_HEX_IRONS_6,
            8 => &IP_HEX_LEGENDRE_8,
            9 => &IP_HEX_WILSON_CORNER_9,
            1_009 => &IP_HEX_WILSON_STABLE_9,
            14 => &IP_HEX_IRONS_14,
            27 => &IP_HEX_LEGENDRE_27,
            _ => return Err("desired number of integration points is not available for Hex class"),
        },
    };
    Ok(ips)
}

//
// Next, we define several constants with integration points data.
//

// -----------------------------------------------------------------------
// -- LIN ----------------------------------------------------------------
// -----------------------------------------------------------------------

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_LIN_LEGENDRE_1: [[f64; 4]; 1] = [
    [0.0, 0.0, 0.0, 2.0]
];

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_LIN_LEGENDRE_2: [[f64; 4]; 2] = [
    [-0.5773502691896257, 0.0, 0.0, 1.0],
    [ 0.5773502691896257, 0.0, 0.0, 1.0],
];

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_LIN_LEGENDRE_3: [[f64; 4]; 3] = [
    [-0.7745966692414834, 0.0, 0.0, 0.5555555555555556],
    [ 0.0000000000000000, 0.0, 0.0, 0.8888888888888888],
    [ 0.7745966692414834, 0.0, 0.0, 0.5555555555555556],
];

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_LIN_LEGENDRE_4: [[f64; 4]; 4] = [
    [-0.8611363115940526, 0.0, 0.0, 0.3478548451374538],
    [-0.3399810435848562, 0.0, 0.0, 0.6521451548625462],
    [ 0.3399810435848562, 0.0, 0.0, 0.6521451548625462],
    [ 0.8611363115940526, 0.0, 0.0, 0.3478548451374538],
];

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_LIN_LEGENDRE_5: [[f64; 4]; 5] = [
    [-0.9061798459386640, 0.0, 0.0, 0.2369268850561891],
    [-0.5384693101056831, 0.0, 0.0, 0.4786286704993665],
    [ 0.0000000000000000, 0.0, 0.0, 0.5688888888888889],
    [ 0.5384693101056831, 0.0, 0.0, 0.4786286704993665],
    [ 0.9061798459386640, 0.0, 0.0, 0.2369268850561891],
];

// -----------------------------------------------------------------------
// -- TRI ----------------------------------------------------------------
// -----------------------------------------------------------------------

/// Internal integration points and weights
#[rustfmt::skip]
pub const IP_TRI_INTERNAL_1: [[f64; 4]; 1] = [
    [1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0 / 2.0],
];

/// Internal integration points and weights
#[rustfmt::skip]
pub const IP_TRI_INTERNAL_3: [[f64; 4]; 3] = [
    [1.0 / 6.0, 1.0 / 6.0, 0.0, 1.0 / 6.0],
    [2.0 / 3.0, 1.0 / 6.0, 0.0, 1.0 / 6.0],
    [1.0 / 6.0, 2.0 / 3.0, 0.0, 1.0 / 6.0],
];

/// Edge integration points and weights
#[rustfmt::skip]
pub const IP_TRI_EDGE_3: [[f64; 4]; 3] = [
    [0.5, 0.5, 0.0, 1.0 / 6.0],
    [0.0, 0.5, 0.0, 1.0 / 6.0],
    [0.5, 0.0, 0.0, 1.0 / 6.0],
];

/// Internal integration points and weights
#[rustfmt::skip]
pub const IP_TRI_INTERNAL_4: [[f64; 4]; 4] = [
    [1.0 / 3.0, 1.0 / 3.0, 0.0, -27.0 / 96.0],
    [1.0 / 5.0, 1.0 / 5.0, 0.0,  25.0 / 96.0],
    [3.0 / 5.0, 1.0 / 5.0, 0.0,  25.0 / 96.0],
    [1.0 / 5.0, 3.0 / 5.0, 0.0,  25.0 / 96.0],
];

/// Internal integration points and weights
#[rustfmt::skip]
pub const IP_TRI_INTERNAL_12: [[f64; 4]; 12] = [
    [0.873821971016996, 0.063089014491502, 0.0, 0.0254224531851035],
    [0.063089014491502, 0.873821971016996, 0.0, 0.0254224531851035],
    [0.063089014491502, 0.063089014491502, 0.0, 0.0254224531851035],
    [0.501426509658179, 0.249286745170910, 0.0, 0.0583931378631895],
    [0.249286745170910, 0.501426509658179, 0.0, 0.0583931378631895],
    [0.249286745170910, 0.249286745170910, 0.0, 0.0583931378631895],
    [0.053145049844817, 0.310352451033784, 0.0, 0.041425537809187 ],
    [0.310352451033784, 0.053145049844817, 0.0, 0.041425537809187 ],
    [0.053145049844817, 0.636502499121398, 0.0, 0.041425537809187 ],
    [0.310352451033784, 0.636502499121398, 0.0, 0.041425537809187 ],
    [0.636502499121398, 0.053145049844817, 0.0, 0.041425537809187 ],
    [0.636502499121398, 0.310352451033784, 0.0, 0.041425537809187 ],
];

/// Internal integration points and weights
#[rustfmt::skip]
pub const IP_TRI_INTERNAL_16: [[f64; 4]; 16] = [
    [3.33333333333333e-01, 3.33333333333333e-01, 0.0, 7.21578038388935e-02],
    [8.14148234145540e-02, 4.59292588292723e-01, 0.0, 4.75458171336425e-02],
    [4.59292588292723e-01, 8.14148234145540e-02, 0.0, 4.75458171336425e-02],
    [4.59292588292723e-01, 4.59292588292723e-01, 0.0, 4.75458171336425e-02],
    [6.58861384496480e-01, 1.70569307751760e-01, 0.0, 5.16086852673590e-02],
    [1.70569307751760e-01, 6.58861384496480e-01, 0.0, 5.16086852673590e-02],
    [1.70569307751760e-01, 1.70569307751760e-01, 0.0, 5.16086852673590e-02],
    [8.98905543365938e-01, 5.05472283170310e-02, 0.0, 1.62292488115990e-02],
    [5.05472283170310e-02, 8.98905543365938e-01, 0.0, 1.62292488115990e-02],
    [5.05472283170310e-02, 5.05472283170310e-02, 0.0, 1.62292488115990e-02],
    [8.39477740995800e-03, 2.63112829634638e-01, 0.0, 1.36151570872175e-02],
    [7.28492392955404e-01, 8.39477740995800e-03, 0.0, 1.36151570872175e-02],
    [2.63112829634638e-01, 7.28492392955404e-01, 0.0, 1.36151570872175e-02],
    [8.39477740995800e-03, 7.28492392955404e-01, 0.0, 1.36151570872175e-02],
    [7.28492392955404e-01, 2.63112829634638e-01, 0.0, 1.36151570872175e-02],
    [2.63112829634638e-01, 8.39477740995800e-03, 0.0, 1.36151570872175e-02],
];

// -----------------------------------------------------------------------
// -- QUA ----------------------------------------------------------------
// -----------------------------------------------------------------------

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_QUA_LEGENDRE_1: [[f64; 4]; 1] = [
    [0.0, 0.0, 0.0, 4.0],
];

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_QUA_LEGENDRE_4: [[f64; 4]; 4] = [
    [-0.5773502691896257, -0.5773502691896257, 0.0, 1.0],
    [ 0.5773502691896257, -0.5773502691896257, 0.0, 1.0],
    [-0.5773502691896257,  0.5773502691896257, 0.0, 1.0],
    [ 0.5773502691896257,  0.5773502691896257, 0.0, 1.0],
];

/// Wilson's integration points and weights. "Corner" version
#[rustfmt::skip]
pub const IP_QUA_WILSON_CORNER_5: [[f64; 4]; 5] = [
    [-1.0, -1.0, 0.0, 0.3333333333333333],
    [ 1.0, -1.0, 0.0, 0.3333333333333333],
    [ 0.0,  0.0, 0.0, 2.6666666666666665],
    [-1.0,  1.0, 0.0, 0.3333333333333333],
    [ 1.0,  1.0, 0.0, 0.3333333333333333],
];

/// Wilson's integration points and weights. "Stable" version version with w0=0.004 and wa=0.999 to mimic 4-point rule
#[rustfmt::skip]
pub const IP_QUA_WILSON_STABLE_5: [[f64; 4]; 5] = [
    [-0.5776391000000000, -0.5776391000000000, 0.0, 0.999],
    [ 0.5776391000000000, -0.5776391000000000, 0.0, 0.999],
    [ 0.0000000000000000,  0.0000000000000000, 0.0, 0.004],
    [-0.5776391000000000,  0.5776391000000000, 0.0, 0.999],
    [ 0.5776391000000000,  0.5776391000000000, 0.0, 0.999],
];

/// Wilson's integration points and weights.
#[rustfmt::skip]
pub const IP_QUA_WILSON_8: [[f64; 4]; 8] = [
    [-0.8819171036881969, -0.8819171036881969, 0.0, 0.1836734693877551],
    [ 0.0000000000000000, -0.6831300510639732, 0.0, 0.8163265306122449],
    [ 0.8819171036881969, -0.8819171036881969, 0.0, 0.1836734693877551],
    [-0.6831300510639732,  0.0000000000000000, 0.0, 0.8163265306122449],
    [ 0.6831300510639732,  0.0000000000000000, 0.0, 0.8163265306122449],
    [-0.8819171036881969,  0.8819171036881969, 0.0, 0.1836734693877551],
    [ 0.0000000000000000,  0.6831300510639732, 0.0, 0.8163265306122449],
    [ 0.8819171036881969,  0.8819171036881969, 0.0, 0.1836734693877551],
];

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_QUA_LEGENDRE_9: [[f64; 4]; 9] = [
    [-0.7745966692414834, -0.7745966692414834, 0.0, 25.0 / 81.0],
    [ 0.0000000000000000, -0.7745966692414834, 0.0, 40.0 / 81.0],
    [ 0.7745966692414834, -0.7745966692414834, 0.0, 25.0 / 81.0],
    [-0.7745966692414834,  0.0000000000000000, 0.0, 40.0 / 81.0],
    [ 0.0000000000000000,  0.0000000000000000, 0.0, 64.0 / 81.0],
    [ 0.7745966692414834,  0.0000000000000000, 0.0, 40.0 / 81.0],
    [-0.7745966692414834,  0.7745966692414834, 0.0, 25.0 / 81.0],
    [ 0.0000000000000000,  0.7745966692414834, 0.0, 40.0 / 81.0],
    [ 0.7745966692414834,  0.7745966692414834, 0.0, 25.0 / 81.0],
];

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_QUA_LEGENDRE_16: [[f64; 4]; 16] = [
    [-0.8611363115940526, -0.8611363115940526, 0.0, 0.1210029932856019],
    [-0.3399810435848563, -0.8611363115940526, 0.0, 0.2268518518518519],
    [ 0.3399810435848563, -0.8611363115940526, 0.0, 0.2268518518518519],
    [ 0.8611363115940526, -0.8611363115940526, 0.0, 0.1210029932856019],
    [-0.8611363115940526, -0.3399810435848563, 0.0, 0.2268518518518519],
    [-0.3399810435848563, -0.3399810435848563, 0.0, 0.4252933030106947],
    [ 0.3399810435848563, -0.3399810435848563, 0.0, 0.4252933030106947],
    [ 0.8611363115940526, -0.3399810435848563, 0.0, 0.2268518518518519],
    [-0.8611363115940526,  0.3399810435848563, 0.0, 0.2268518518518519],
    [-0.3399810435848563,  0.3399810435848563, 0.0, 0.4252933030106947],
    [ 0.3399810435848563,  0.3399810435848563, 0.0, 0.4252933030106947],
    [ 0.8611363115940526,  0.3399810435848563, 0.0, 0.2268518518518519],
    [-0.8611363115940526,  0.8611363115940526, 0.0, 0.1210029932856019],
    [-0.3399810435848563,  0.8611363115940526, 0.0, 0.2268518518518519],
    [ 0.3399810435848563,  0.8611363115940526, 0.0, 0.2268518518518519],
    [ 0.8611363115940526,  0.8611363115940526, 0.0, 0.1210029932856019],
];

// -----------------------------------------------------------------------
// -- TET ----------------------------------------------------------------
// -----------------------------------------------------------------------

/// Internal integration points and weights, 1 point, degree 1
/// 
/// Reference: Wriggers (2008), page 12, table 4.5,
/// Nonlinear Finite Element Methods, Springer.
#[rustfmt::skip]
pub const IP_TET_INTERNAL_1: [[f64; 4]; 1] = [
    [1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 6.0],
];

/// Internal integration points and weights, 4 points, degree 2
/// 
/// References: Felippa Advanced FEM, page 17-21, Table 17.1. See also Felippa (2004) A compendium of
/// FEM integration formulas for symbolic work Engineering Computations, 21, 7/8, pg 867
///
/// See Mathematica code at the bottom of this source file and also in the directory `data/derivations`
#[rustfmt::skip]
pub const IP_TET_INTERNAL_4: [[f64; 4]; 4] = [
    [0.13819660112501051518, 0.13819660112501051518, 0.13819660112501051518, 0.041666666666666666667],
    [0.58541019662496845446, 0.13819660112501051518, 0.13819660112501051518, 0.041666666666666666667], 
    [0.13819660112501051518, 0.58541019662496845446, 0.13819660112501051518, 0.041666666666666666667],
    [0.13819660112501051518, 0.13819660112501051518, 0.58541019662496845446, 0.041666666666666666667],
];

/// Internal integration points and weights, 4 points, degree 3
/// 
/// Reference: Wriggers (2008), page 12, table 4.5,
/// Nonlinear Finite Element Methods, Springer.
#[rustfmt::skip]
pub const IP_TET_INTERNAL_5: [[f64; 4]; 5] = [
    [1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, -2.0 / 15.0],
    [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0,  3.0 / 40.0],
    [1.0 / 6.0, 1.0 / 6.0, 1.0 / 2.0,  3.0 / 40.0],
    [1.0 / 6.0, 1.0 / 2.0, 1.0 / 6.0,  3.0 / 40.0],
    [1.0 / 2.0, 1.0 / 6.0, 1.0 / 6.0,  3.0 / 40.0],
];

/// Internal integration points and weights, 8 points, degree 3
/// 
/// References: Felippa Advanced FEM, page 17-21, Table 17.1. See also Felippa (2004) A compendium of
/// FEM integration formulas for symbolic work Engineering Computations, 21, 7/8, pg 867
///
/// See Mathematica code at the bottom of this source file and also in the directory `data/derivations`
#[rustfmt::skip]
pub const IP_TET_FELIPPA_8: [[f64; 4]; 8] = [
    [0.32805469671142664734,  0.32805469671142664734,  0.32805469671142664734,  0.023087994418643690387], 
    [0.015835909865720057993, 0.32805469671142664734,  0.32805469671142664734,  0.023087994418643690387], 
    [0.32805469671142664734,  0.015835909865720057993, 0.32805469671142664734,  0.023087994418643690387], 
    [0.32805469671142664734,  0.32805469671142664734,  0.015835909865720057993, 0.023087994418643690387], 
    [0.10695227393293068277,  0.10695227393293068277,  0.10695227393293068277,  0.018578672248022976279], 
    [0.67914317820120795168,  0.10695227393293068277,  0.10695227393293068277,  0.018578672248022976279], 
    [0.10695227393293068277,  0.67914317820120795168,  0.10695227393293068277,  0.018578672248022976279], 
    [0.10695227393293068277,  0.10695227393293068277,  0.67914317820120795168,  0.018578672248022976279],
];

/// Internal integration points and weights, 14 points, degree 4
/// 
/// References: Felippa Advanced FEM, page 17-21, Table 17.1. See also Felippa (2004) A compendium of
/// FEM integration formulas for symbolic work Engineering Computations, 21, 7/8, pg 867
///
/// See Mathematica code at the bottom of this source file and also in the directory `data/derivations`
#[rustfmt::skip]
pub const IP_TET_FELIPPA_14: [[f64; 4]; 14] = [
    [0.092735250310891226402, 0.092735250310891226402, 0.092735250310891226402, 0.012248840519393658257 ], 
    [0.72179424906732632079,  0.092735250310891226402, 0.092735250310891226402, 0.012248840519393658257 ], 
    [0.092735250310891226402, 0.72179424906732632079,  0.092735250310891226402, 0.012248840519393658257 ], 
    [0.092735250310891226402, 0.092735250310891226402, 0.72179424906732632079,  0.012248840519393658257 ], 
    [0.31088591926330060980,  0.31088591926330060980,  0.31088591926330060980,  0.018781320953002641800 ], 
    [0.067342242210098170608, 0.31088591926330060980,  0.31088591926330060980,  0.018781320953002641800 ], 
    [0.31088591926330060980,  0.067342242210098170608, 0.31088591926330060980,  0.018781320953002641800 ], 
    [0.31088591926330060980,  0.31088591926330060980,  0.067342242210098170608, 0.018781320953002641800 ], 
    [0.045503704125649649492, 0.45449629587435035051,  0.45449629587435035051,  0.0070910034628469110730], 
    [0.45449629587435035051,  0.045503704125649649492, 0.45449629587435035051,  0.0070910034628469110730], 
    [0.45449629587435035051,  0.45449629587435035051,  0.045503704125649649492, 0.0070910034628469110730], 
    [0.045503704125649649492, 0.045503704125649649492, 0.45449629587435035051,  0.0070910034628469110730], 
    [0.045503704125649649492, 0.45449629587435035051,  0.045503704125649649492, 0.0070910034628469110730], 
    [0.45449629587435035051,  0.045503704125649649492, 0.045503704125649649492, 0.0070910034628469110730],
];

/// Internal integration points and weights, 15 points, degree 5
/// 
/// References: Felippa Advanced FEM, page 17-21, Table 17.1. See also Felippa (2004) A compendium of
/// FEM integration formulas for symbolic work Engineering Computations, 21, 7/8, pg 867
///
/// See Mathematica code at the bottom of this source file and also in the directory `data/derivations`
#[rustfmt::skip]
pub const IP_TET_FELIPPA_15: [[f64; 4]; 15] = [
    [0.091971078052723032789, 0.091971078052723032789, 0.091971078052723032789, 0.011989513963169770002 ], 
    [0.72408676584183090163,  0.091971078052723032789, 0.091971078052723032789, 0.011989513963169770002 ], 
    [0.091971078052723032789, 0.72408676584183090163,  0.091971078052723032789, 0.011989513963169770002 ], 
    [0.091971078052723032789, 0.091971078052723032789, 0.72408676584183090163,  0.011989513963169770002 ], 
    [0.31979362782962990839,  0.31979362782962990839,  0.31979362782962990839,  0.011511367871045397547 ], 
    [0.040619116511110274837, 0.31979362782962990839,  0.31979362782962990839,  0.011511367871045397547 ], 
    [0.31979362782962990839,  0.040619116511110274837, 0.31979362782962990839,  0.011511367871045397547 ], 
    [0.31979362782962990839,  0.31979362782962990839,  0.040619116511110274837, 0.011511367871045397547 ], 
    [0.44364916731037084426,  0.056350832689629155741, 0.056350832689629155741, 0.0088183421516754850088], 
    [0.056350832689629155741, 0.44364916731037084426,  0.056350832689629155741, 0.0088183421516754850088], 
    [0.056350832689629155741, 0.056350832689629155741, 0.44364916731037084426,  0.0088183421516754850088], 
    [0.44364916731037084426,  0.44364916731037084426,  0.056350832689629155741, 0.0088183421516754850088], 
    [0.44364916731037084426,  0.056350832689629155741, 0.44364916731037084426,  0.0088183421516754850088], 
    [0.056350832689629155741, 0.44364916731037084426,  0.44364916731037084426,  0.0088183421516754850088], 
    [0.25000000000000000000,  0.25000000000000000000,  0.25000000000000000000,  0.019753086419753086420 ],
];

/// Internal integration points and weights, 24 points, degree 6
/// 
/// References: Felippa Advanced FEM, page 17-21, Table 17.1. See also Felippa (2004) A compendium of
/// FEM integration formulas for symbolic work Engineering Computations, 21, 7/8, pg 867
///
/// See Mathematica code at the bottom of this source file and also in the directory `data/derivations`
#[rustfmt::skip]
pub const IP_TET_FELIPPA_24: [[f64; 4]; 24] = [
    [0.21460287125915202929,  0.21460287125915202929,  0.21460287125915202929,  0.0066537917096945820166], 
    [0.35619138622254391213,  0.21460287125915202929,  0.21460287125915202929,  0.0066537917096945820166], 
    [0.21460287125915202929,  0.35619138622254391213,  0.21460287125915202929,  0.0066537917096945820166], 
    [0.21460287125915202929,  0.21460287125915202929,  0.35619138622254391213,  0.0066537917096945820166], 
    [0.040673958534611353116, 0.040673958534611353116, 0.040673958534611353116, 0.0016795351758867738247], 
    [0.87797812439616594065,  0.040673958534611353116, 0.040673958534611353116, 0.0016795351758867738247], 
    [0.040673958534611353116, 0.87797812439616594065,  0.040673958534611353116, 0.0016795351758867738247], 
    [0.040673958534611353116, 0.040673958534611353116, 0.87797812439616594065,  0.0016795351758867738247], 
    [0.32233789014227551034,  0.32233789014227551034,  0.32233789014227551034,  0.0092261969239424536825], 
    [0.032986329573173468968, 0.32233789014227551034,  0.32233789014227551034,  0.0092261969239424536825], 
    [0.32233789014227551034,  0.032986329573173468968, 0.32233789014227551034,  0.0092261969239424536825], 
    [0.32233789014227551034,  0.32233789014227551034,  0.032986329573173468968, 0.0092261969239424536825], 
    [0.26967233145831580803,  0.063661001875017525299, 0.063661001875017525299, 0.0080357142857142857143], 
    [0.063661001875017525299, 0.26967233145831580803,  0.063661001875017525299, 0.0080357142857142857143], 
    [0.063661001875017525299, 0.063661001875017525299, 0.26967233145831580803,  0.0080357142857142857143], 
    [0.60300566479164914137,  0.26967233145831580803,  0.063661001875017525299, 0.0080357142857142857143], 
    [0.60300566479164914137,  0.063661001875017525299, 0.26967233145831580803,  0.0080357142857142857143], 
    [0.063661001875017525299, 0.60300566479164914137,  0.26967233145831580803,  0.0080357142857142857143], 
    [0.60300566479164914137,  0.063661001875017525299, 0.063661001875017525299, 0.0080357142857142857143], 
    [0.063661001875017525299, 0.60300566479164914137,  0.063661001875017525299, 0.0080357142857142857143], 
    [0.063661001875017525299, 0.063661001875017525299, 0.60300566479164914137,  0.0080357142857142857143], 
    [0.26967233145831580803,  0.60300566479164914137,  0.063661001875017525299, 0.0080357142857142857143], 
    [0.26967233145831580803,  0.063661001875017525299, 0.60300566479164914137,  0.0080357142857142857143], 
    [0.063661001875017525299, 0.26967233145831580803,  0.60300566479164914137,  0.0080357142857142857143],
];

// -----------------------------------------------------------------------
// -- HEX ----------------------------------------------------------------
// -----------------------------------------------------------------------

/// Iron's integration points and weights
#[rustfmt::skip]
pub const IP_HEX_IRONS_6: [[f64; 4]; 6] = [
    [-1.0,  0.0,  0.0, 4.0 / 3.0],
    [ 1.0,  0.0,  0.0, 4.0 / 3.0],
    [ 0.0, -1.0,  0.0, 4.0 / 3.0],
    [ 0.0,  1.0,  0.0, 4.0 / 3.0],
    [ 0.0,  0.0, -1.0, 4.0 / 3.0],
    [ 0.0,  0.0,  1.0, 4.0 / 3.0],
];

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_HEX_LEGENDRE_8: [[f64; 4]; 8] = [
    [-0.5773502691896257, -0.5773502691896257, -0.5773502691896257, 1.0],
    [ 0.5773502691896257, -0.5773502691896257, -0.5773502691896257, 1.0],
    [-0.5773502691896257,  0.5773502691896257, -0.5773502691896257, 1.0],
    [ 0.5773502691896257,  0.5773502691896257, -0.5773502691896257, 1.0],
    [-0.5773502691896257, -0.5773502691896257,  0.5773502691896257, 1.0],
    [ 0.5773502691896257, -0.5773502691896257,  0.5773502691896257, 1.0],
    [-0.5773502691896257,  0.5773502691896257,  0.5773502691896257, 1.0],
    [ 0.5773502691896257,  0.5773502691896257,  0.5773502691896257, 1.0],
];

/// Wilson's integration points and weights. "Corner" version
#[rustfmt::skip]
pub const IP_HEX_WILSON_CORNER_9: [[f64; 4]; 9] = [
    [-1.0, -1.0, -1.0,  0.3333333333333333],
    [ 1.0, -1.0, -1.0,  0.3333333333333333],
    [-1.0,  1.0, -1.0,  0.3333333333333333],
    [ 1.0,  1.0, -1.0,  0.3333333333333333],
    [ 0.0,  0.0,  0.0,  5.3333333333333330],
    [-1.0, -1.0,  1.0,  0.3333333333333333],
    [ 1.0, -1.0,  1.0,  0.3333333333333333],
    [-1.0,  1.0,  1.0,  0.3333333333333333],
    [ 1.0,  1.0,  1.0,  0.3333333333333333],
];

/// Wilson's integration points and weights. "Stable" version
#[rustfmt::skip]
pub const IP_HEX_WILSON_STABLE_9: [[f64; 4]; 9] = [
    [-0.5776391, -0.5776391, -0.5776391000000000,  0.999],
    [ 0.5776391, -0.5776391, -0.5776391000000000,  0.999],
    [-0.5776391,  0.5776391, -0.5776391000000000,  0.999],
    [ 0.5776391,  0.5776391, -0.5776391000000000,  0.999],
    [ 0.0000000,  0.0000000,  0.0000000000000000,  0.008],
    [-0.5776391, -0.5776391,  0.5776391000000000,  0.999],
    [ 0.5776391, -0.5776391,  0.5776391000000000,  0.999],
    [-0.5776391,  0.5776391,  0.5776391000000000,  0.999],
    [ 0.5776391,  0.5776391,  0.5776391000000000,  0.999],
];

/// Iron's integration points and weights
#[rustfmt::skip]
pub const IP_HEX_IRONS_14: [[f64; 4]; 14] = [
    [ 0.7958224257542215,  0.0000000000000000,  0.0000000000000000,  0.8864265927977839],
    [-0.7958224257542215,  0.0000000000000000,  0.0000000000000000,  0.8864265927977839],
    [ 0.0000000000000000,  0.7958224257542215,  0.0000000000000000,  0.8864265927977839],
    [ 0.0000000000000000, -0.7958224257542215,  0.0000000000000000,  0.8864265927977839],
    [ 0.0000000000000000,  0.0000000000000000,  0.7958224257542215,  0.8864265927977839],
    [ 0.0000000000000000,  0.0000000000000000, -0.7958224257542215,  0.8864265927977839],
    [ 0.7587869106393281,  0.7587869106393281,  0.7587869106393281,  0.3351800554016621],
    [-0.7587869106393281,  0.7587869106393281,  0.7587869106393281,  0.3351800554016621],
    [ 0.7587869106393281, -0.7587869106393281,  0.7587869106393281,  0.3351800554016621],
    [-0.7587869106393281, -0.7587869106393281,  0.7587869106393281,  0.3351800554016621],
    [ 0.7587869106393281,  0.7587869106393281, -0.7587869106393281,  0.3351800554016621],
    [-0.7587869106393281,  0.7587869106393281, -0.7587869106393281,  0.3351800554016621],
    [ 0.7587869106393281, -0.7587869106393281, -0.7587869106393281,  0.3351800554016621],
    [-0.7587869106393281, -0.7587869106393281, -0.7587869106393281,  0.3351800554016621],
];

/// Conventional Legendre integration points and weights
#[rustfmt::skip]
pub const IP_HEX_LEGENDRE_27: [[f64; 4]; 27] = [
    [-0.774596669241483, -0.774596669241483, -0.774596669241483, 0.171467764060357],
    [ 0.000000000000000, -0.774596669241483, -0.774596669241483, 0.274348422496571],
    [ 0.774596669241483, -0.774596669241483, -0.774596669241483, 0.171467764060357],
    [-0.774596669241483,  0.000000000000000, -0.774596669241483, 0.274348422496571],
    [ 0.000000000000000,  0.000000000000000, -0.774596669241483, 0.438957475994513],
    [ 0.774596669241483,  0.000000000000000, -0.774596669241483, 0.274348422496571],
    [-0.774596669241483,  0.774596669241483, -0.774596669241483, 0.171467764060357],
    [ 0.000000000000000,  0.774596669241483, -0.774596669241483, 0.274348422496571],
    [ 0.774596669241483,  0.774596669241483, -0.774596669241483, 0.171467764060357],
    [-0.774596669241483, -0.774596669241483,  0.000000000000000, 0.274348422496571],
    [ 0.000000000000000, -0.774596669241483,  0.000000000000000, 0.438957475994513],
    [ 0.774596669241483, -0.774596669241483,  0.000000000000000, 0.274348422496571],
    [-0.774596669241483,  0.000000000000000,  0.000000000000000, 0.438957475994513],
    [ 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.702331961591221],
    [ 0.774596669241483,  0.000000000000000,  0.000000000000000, 0.438957475994513],
    [-0.774596669241483,  0.774596669241483,  0.000000000000000, 0.274348422496571],
    [ 0.000000000000000,  0.774596669241483,  0.000000000000000, 0.438957475994513],
    [ 0.774596669241483,  0.774596669241483,  0.000000000000000, 0.274348422496571],
    [-0.774596669241483, -0.774596669241483,  0.774596669241483, 0.171467764060357],
    [ 0.000000000000000, -0.774596669241483,  0.774596669241483, 0.274348422496571],
    [ 0.774596669241483, -0.774596669241483,  0.774596669241483, 0.171467764060357],
    [-0.774596669241483,  0.000000000000000,  0.774596669241483, 0.274348422496571],
    [ 0.000000000000000,  0.000000000000000,  0.774596669241483, 0.438957475994513],
    [ 0.774596669241483,  0.000000000000000,  0.774596669241483, 0.274348422496571],
    [-0.774596669241483,  0.774596669241483,  0.774596669241483, 0.171467764060357],
    [ 0.000000000000000,  0.774596669241483,  0.774596669241483, 0.274348422496571],
    [ 0.774596669241483,  0.774596669241483,  0.774596669241483, 0.171467764060357],
];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{default_integ_points, select_integ_points};
    use crate::shapes::{GeoClass, GeoKind};
    use crate::StrError;

    #[test]
    fn default_integ_points_works() {
        // Lin
        assert_eq!(default_integ_points(GeoKind::Lin2).len(), 2);
        assert_eq!(default_integ_points(GeoKind::Lin3).len(), 3);
        assert_eq!(default_integ_points(GeoKind::Lin4).len(), 4);
        assert_eq!(default_integ_points(GeoKind::Lin5).len(), 5);
        // Tri
        assert_eq!(default_integ_points(GeoKind::Tri3).len(), 1);
        assert_eq!(default_integ_points(GeoKind::Tri6).len(), 4);
        assert_eq!(default_integ_points(GeoKind::Tri10).len(), 12);
        assert_eq!(default_integ_points(GeoKind::Tri15).len(), 16);
        // Qua
        assert_eq!(default_integ_points(GeoKind::Qua4).len(), 4);
        assert_eq!(default_integ_points(GeoKind::Qua8).len(), 9);
        assert_eq!(default_integ_points(GeoKind::Qua9).len(), 9);
        assert_eq!(default_integ_points(GeoKind::Qua12).len(), 16);
        assert_eq!(default_integ_points(GeoKind::Qua16).len(), 16);
        assert_eq!(default_integ_points(GeoKind::Qua17).len(), 16);
        // Tet
        assert_eq!(default_integ_points(GeoKind::Tet4).len(), 1);
        assert_eq!(default_integ_points(GeoKind::Tet10).len(), 14);
        // Hex
        assert_eq!(default_integ_points(GeoKind::Hex8).len(), 8);
        assert_eq!(default_integ_points(GeoKind::Hex20).len(), 27);
    }

    #[test]
    fn select_integ_points_works() -> Result<(), StrError> {
        // Lin
        for n_integ_point in [1, 2, 3, 4, 5] {
            let ips = select_integ_points(GeoClass::Lin, n_integ_point)?;
            assert_eq!(ips.len(), n_integ_point);
        }
        assert_eq!(
            select_integ_points(GeoClass::Lin, 100).err(),
            Some("desired number of integration points is not available for Lin class")
        );

        // Tri
        for n_integ_point in [1, 3, 4, 12, 16] {
            let ips = select_integ_points(GeoClass::Tri, n_integ_point)?;
            assert_eq!(ips.len(), n_integ_point);
        }
        let ips = select_integ_points(GeoClass::Tri, 1_003)?;
        assert_eq!(ips.len(), 3);
        assert_eq!(
            select_integ_points(GeoClass::Tri, 100).err(),
            Some("desired number of integration points is not available for Tri class")
        );

        // Qua
        for n_integ_point in [1, 4, 5, 8, 9, 16] {
            let ips = select_integ_points(GeoClass::Qua, n_integ_point)?;
            assert_eq!(ips.len(), n_integ_point);
        }
        let ips = select_integ_points(GeoClass::Qua, 1_005)?;
        assert_eq!(ips.len(), 5);
        assert_eq!(
            select_integ_points(GeoClass::Qua, 100).err(),
            Some("desired number of integration points is not available for Qua class")
        );

        // Tet
        for n_integ_point in [1, 4, 5, 8, 14] {
            let ips = select_integ_points(GeoClass::Tet, n_integ_point)?;
            assert_eq!(ips.len(), n_integ_point);
        }
        assert_eq!(
            select_integ_points(GeoClass::Tet, 100).err(),
            Some("desired number of integration points is not available for Tet class")
        );

        // Hex
        for n_integ_point in [6, 8, 9, 14, 27] {
            let ips = select_integ_points(GeoClass::Hex, n_integ_point)?;
            assert_eq!(ips.len(), n_integ_point);
        }
        let ips = select_integ_points(GeoClass::Hex, 1_009)?;
        assert_eq!(ips.len(), 9);
        assert_eq!(
            select_integ_points(GeoClass::Hex, 100).err(),
            Some("desired number of integration points is not available for Hex class")
        );
        Ok(())
    }
}

// Modified Felippa's Mathematica Code
//
// From Advanced Finite Elements, Chapter 17, page 17-22
//
// See also directory data/derivations
//
// The modification is just a division by 6 near the last line.
// In this way, the weights are consistent with other IPs here.
//
// ```mathematica
// TetrGaussRuleInfo[{rule_, numer_}, point_] :=
//   Module[{jk6 = {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}},
//     jk12 = {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}, {2,
//        1}, {3, 1}, {4, 1}, {3, 2}, {4, 2}, {4, 3}}, i = point, j, k,
//     g1, g2, g3, g4, h1, w1, w2, w3, eps = 10.^(-16),
//     info = {{Null, Null, Null, Null}, 0}},
//    If[rule == 1, info = {{1/4, 1/4, 1/4, 1/4}, 1}];
//    If[rule == 4, g1 = (5 - Sqrt[5])/20; h1 = (5 + 3*Sqrt[5])/20;
//     info = {{g1, g1, g1, g1}, 1/4}; info[[1, i]] = h1];
//    If[rule == 8, j = i - 4;
//     g1 = (55 - 3*Sqrt[17] + Sqrt[1022 - 134*Sqrt[17]])/196;
//     g2 = (55 - 3*Sqrt[17] - Sqrt[1022 - 134*Sqrt[17]])/196;
//     w1 = 1/8 + Sqrt[(1715161837 - 406006699*Sqrt[17])/23101]/3120;
//     w2 = 1/8 - Sqrt[(1715161837 - 406006699*Sqrt[17])/23101]/3120;
//     If[j <= 0, info = {{g1, g1, g1, g1}, w1}; info[[1, i]] = 1 - 3*g1];
//     If[j > 0, info = {{g2, g2, g2, g2}, w2}; info[[1, j]] = 1 - 3*g2]];
//    If[rule == -8, j = i - 4;
//     If[j <= 0, info = {{0, 0, 0, 0}, 1/40}; info[[1, i]] = 1];
//     If[j > 0, info = {{1, 1, 1, 1}/3, 9/40}; info[[1, j]] = 0]];
//    If[rule == 14,(*g1,g2+roots of P(g)=0,P=9+96*g-1712*g^2-30464*
//     g^3-127232*g^4+86016*g^5+1060864*g^6*)
//     g1 = 0.09273525031089122640232391373703060;
//     g2 = 0.31088591926330060979734573376345783;
//     g3 = 0.45449629587435035050811947372066056;
//     If[! numer, {g1, g2, g3} = Rationalize[{g1, g2, g3}, eps]];
//     w1 = (-1 + 6*g2*(2 + g2*(-7 + 8*g2)) + 14*g3 -
//         60*g2*(3 + 4*g2*(-3 + 4*g2))*g3 +
//         4*(-7 + 30*g2*(3 + 4*g2*(-3 + 4*g2)))*
//          g3^2)/(120*(g1 - g2)*(g2*(-3 + 8*g2) + 6*g3 +
//           8*g2*(-3 + 4*g2)*g3 - 4*(3 + 4*g2*(-3 + 4*g2))*g3^2 +
//           8*g1^2*(1 + 12*g2*(-1 + 2*g2) + 4*g3 - 8*g3^2) +
//           g1*(-3 - 96*g2^2 + 24*g3*(-1 + 2*g3) +
//              g2*(44 + 32*(1 - 2*g3)*g3))));
//     w2 = (-1 - 20*(1 + 12*g1*(2*g1 - 1))*w1 +
//         20*g3*(2*g3 - 1)*(4*w1 - 1))/(20*(1 + 12*g2*(2*g2 - 1) +
//           4*g3 - 8*g3^2));
//     If[i < 5, info = {{g1, g1, g1, g1}, w1};
//      info[[1, i]] = 1 - 3*g1];
//     If[i > 4 && i < 9, info = {{g2, g2, g2, g2}, w2};
//      info[[1, i - 4]] = 1 - 3*g2];
//     If[i > 8, info = {{g3, g3, g3, g3}, 1/6 - 2*(w1 + w2)/3};
//      {j, k} = jk6[[i - 8]]; info[[1, j]] = info[[1, k]] = 1/2 - g3]];
//    If[rule == -14,
//     g1 = (243 - 51*Sqrt[11] + 2*Sqrt[16486 - 9723*Sqrt[11]/2])/356;
//     g2 = (243 - 51*Sqrt[11] - 2*Sqrt[16486 - 9723*Sqrt[11]/2])/356;
//     w1 = 31/280 + Sqrt[(13686301 - 3809646*Sqrt[11])/5965]/600;
//     w2 = 31/280 - Sqrt[(13686301 - 3809646*Sqrt[11])/5965]/600;
//     If[i < 5, info = {{g1, g1, g1, g1}, w1};
//      info[[1, i]] = 1 - 3*g1];
//     If[i > 4 && i < 9, info = {{g2, g2, g2, g2}, w2};
//      info[[1, i - 4]] = 1 - 3*g2];
//     If[i > 8 && i < 15, info = {{0, 0, 0, 0}, 2/105};
//      {j, k} = jk6[[i - 8]]; info[[1, j]] = info[[1, k]] = 1/2]];
//    If[rule == 15, g1 = (7 - Sqrt[15])/34; g2 = 7/17 - g1;
//     g3 = (10 - 2*Sqrt[15])/40;
//     w1 = (2665 + 14*Sqrt[15])/37800; w2 = (2665 - 14*Sqrt[15])/37800;
//     If[i < 5, info = {{g1, g1, g1, g1}, w1};
//      info[[1, i]] = 1 - 3*g1];
//     If[i > 4 && i < 9, info = {{g2, g2, g2, g2}, w2};
//      info[[1, i - 4]] = 1 - 3*g2];
//     If[i > 8 && i < 15, info = {{g3, g3, g3, g3}, 10/189};
//      {j, k} = jk6[[i - 8]]; info[[1, j]] = info[[1, k]] = 1/2 - g3];
//     If[i == 15, info = {{1/4, 1/4, 1/4, 1/4}, 16/135}]];
//    If[rule == -15, g1 = (13 - Sqrt[91])/52;
//     If[i < 5, info = {{1, 1, 1, 1}/3, 81/2240}; info[[1, i]] = 0];
//     If[i > 4 && i < 9, info = {{1, 1, 1, 1}/11, 161051/2304960};
//      info[[1, i - 4]] = 8/11];
//     If[i > 8 && i < 15, info = {{g1, g1, g1, g1}, 338/5145};
//      {j, k} = jk6[[i - 8]]; info[[1, j]] = info[[1, k]] = 1/2 - g1];
//     If[i == 15, info = {{1/4, 1/4, 1/4, 1/4}, 6544/36015}]];
//    If[rule == 24, g1 = 0.214602871259152029288839219386284991;
//     g2 = 0.040673958534611353115579448956410059;
//     g3 = 0.322337890142275510343994470762492125;
//     If[! numer, {g1, g2, g3} = Rationalize[{g1, g2, g3}, eps]];
//     w1 = (85 + 2*g2*(-319 + 9*Sqrt[5] + 624*g2) - 638*g3 -
//         24*g2*(-229 + 472*g2)*g3 +
//         96*(13 + 118*g2*(-1 + 2*g2))*g3^2 +
//         9*Sqrt[5]*(-1 + 2*g3))/(13440*(g1 - g2)*(g1 - g3)*(3 - 8*g2 +
//           8*g1*(-1 + 2*g2) - 8*g3 + 16*(g1 + g2)*g3));
//     w2 = -(85 + 2*g1*(-319 + 9*Sqrt[5] + 624*g1) - 638*g3 -
//          24*g1*(-229 + 472*g1)*g3 +
//          96*(13 + 118*g1*(-1 + 2*g1))*g3^2 +
//          9*Sqrt[5]*(-1 + 2*g3))/(13440*(g1 - g2)*(g2 - g3)*(3 -
//           8*g2 + 8*g1*(-1 + 2*g2) - 8*g3 + 16*(g1 + g2)*g3));
//     w3 = (85 + 2*g1*(-319 + 9*Sqrt[5] + 624*g1) - 638*g2 -
//         24*g1*(-229 + 472*g1)*g2 +
//         96*(13 + 118*g1*(-1 + 2*g1))*g2^2 +
//         9*Sqrt[5]*(-1 + 2*g2))/(13440*(g1 - g3)*(g2 - g3)*(3 - 8*g2 +
//           8*g1*(-1 + 2*g2) - 8*g3 + 16*(g1 + g2)*g3));
//     g4 = (3 - Sqrt[5])/12; h4 = (5 + Sqrt[5])/12;
//     p4 = (1 + Sqrt[5])/12;
//     If[i < 5, info = {{g1, g1, g1, g1}, w1};
//      info[[1, i]] = 1 - 3*g1];
//     If[i > 4 && i < 9, info = {{g2, g2, g2, g2}, w2};
//      info[[1, i - 4]] = 1 - 3*g2];
//     If[i > 8 && i < 13, info = {{g3, g3, g3, g3}, w3};
//      info[[1, i - 8]] = 1 - 3*g3];
//     If[i > 12, info = {{g4, g4, g4, g4}, 27/560};
//      {j, k} = jk12[[i - 12]]; info[[1, j]] = h4; info[[1, k]] = p4]];
//    info[[2]] = info[[2]]/6;(*we include the division by 6 directly here *)
//    If[numer, Return[N[info, 20]], Return[Simplify[info]]];
// ];
// ```

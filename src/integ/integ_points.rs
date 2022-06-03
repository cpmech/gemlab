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
/// * `4` -- Internal integration points and weights, degree 2 (from Felippa's code)
/// * `5` -- Internal integration points and weights, degree 3 (with negative weight)
/// * `8` -- Internal (Felippa's) ips and weights, degree 3 (from Felippa's code)
/// * `14` -- Internal (Felippa's) ips and weights, degree 14 (from Felippa's code)
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
/// Reference: Felippa Advanced FEM, page 17-21, Table 17.1
/// See also Felippa (2004) A compendium of FEM integration formulas for symbolic work
/// Engineering Computations, 21, 7/8, pg 867
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
/// Reference: Felippa Advanced FEM, page 17-21, Table 17.1
/// See also Felippa (2004) A compendium of FEM integration formulas for symbolic work
/// Engineering Computations, 21, 7/8, pg 867
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
/// References: Felippa Advanced FEM, page 17-21, Table 17.1
/// See also Felippa (2004) A compendium of FEM integration formulas for symbolic work
/// Engineering Computations, 21, 7/8, pg 867
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

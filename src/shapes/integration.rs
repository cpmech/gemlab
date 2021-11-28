use super::*;
use crate::StrError;

/// Defines options to select integration points
#[derive(Eq, PartialEq)]
pub enum OptionIntPoint {
    /// Option for Tri (points on edge)
    Edge,

    /// Option for Qua and Hex (Wilson formula)
    WilsonStable,
}

impl Shape {
    /// Selects integrations points and weights
    ///
    /// # Options
    ///
    /// ## nip for Lin class
    ///
    /// * `1` -- Legendre points
    /// * `2` -- Legendre points
    /// * `3` -- Legendre points
    /// * `4` -- Legendre points
    /// * `5` -- Legendre points
    ///
    /// ## nip for Tri class
    ///
    /// * `1` -- Internal points
    /// * `3`:
    ///     - `Edge` -- Points on edge
    ///     - Otherwise, internal points
    /// * `4` -- Internal points
    /// * `12` -- Internal points
    /// * `16` -- Internal points
    ///
    /// ## nip for Qua class
    ///
    /// * `1` -- Legendre points
    /// * `4` -- Legendre points
    /// * `5`:
    ///     - `WilsonStable` -- Wilson's "Stable" version with w0=0.004 and wa=0.999 to mimic 4-point rule
    ///     - Otherwise, Wilson's formula
    /// * `9` -- Legendre points
    /// * `16` -- Legendre points
    ///
    /// ## nip for Tet class
    ///
    /// * `1` -- Internal points
    /// * `4` -- Internal points
    /// * `5` -- Internal points
    /// * `6` -- Internal points
    ///
    /// ## nip for Hex class
    ///
    /// * `6` -- Iron's formula
    /// * `8` -- Legendre points
    /// * `9`:
    ///     - `WilsonStable` -- Wilson's "Stable" version
    ///     - Otherwise, Wilson's formula
    /// * `14` -- Iron's formula
    /// * `27` -- Legendre points
    pub fn select_int_points(&mut self, nip: usize, option: OptionIntPoint) -> Result<(), StrError> {
        match self.class {
            GeoClass::Lin => match nip {
                1 => self.ip_data = &IP_LIN_LEGENDRE_1,
                2 => self.ip_data = &IP_LIN_LEGENDRE_2,
                3 => self.ip_data = &IP_LIN_LEGENDRE_3,
                4 => self.ip_data = &IP_LIN_LEGENDRE_4,
                5 => self.ip_data = &IP_LIN_LEGENDRE_5,
                _ => return Err("number of integration points is not available for Lin class"),
            },
            GeoClass::Tri => match nip {
                1 => self.ip_data = &IP_TRI_INTERNAL_1,
                3 => {
                    if option == OptionIntPoint::Edge {
                        self.ip_data = &IP_TRI_EDGE_3
                    } else {
                        self.ip_data = &IP_TRI_INTERNAL_3
                    }
                }
                4 => self.ip_data = &IP_TRI_INTERNAL_4,
                12 => self.ip_data = &IP_TRI_INTERNAL_12,
                16 => self.ip_data = &IP_TRI_INTERNAL_16,
                _ => return Err("number of integration points is not available for Tri class"),
            },
            GeoClass::Qua => match nip {
                1 => self.ip_data = &IP_QUA_LEGENDRE_1,
                4 => self.ip_data = &IP_QUA_LEGENDRE_4,
                5 => {
                    if option == OptionIntPoint::WilsonStable {
                        self.ip_data = &IP_QUA_WILSON_STABLE_5
                    } else {
                        self.ip_data = &IP_QUA_WILSON_CORNER_5
                    }
                }
                8 => self.ip_data = &IP_QUA_WILSON_8,
                9 => self.ip_data = &IP_QUA_LEGENDRE_9,
                16 => self.ip_data = &IP_QUA_LEGENDRE_16,
                _ => return Err("number of integration points is not available for Qua class"),
            },
            GeoClass::Tet => match nip {
                1 => self.ip_data = &IP_TET_INTERNAL_1,
                4 => self.ip_data = &IP_TET_INTERNAL_4,
                5 => self.ip_data = &IP_TET_INTERNAL_5,
                6 => self.ip_data = &IP_TET_INTERNAL_6,
                _ => return Err("number of integration points is not available for Tet class"),
            },
            GeoClass::Hex => match nip {
                6 => self.ip_data = &IP_HEX_IRONS_6,
                8 => self.ip_data = &IP_HEX_LEGENDRE_8,
                9 => {
                    if option == OptionIntPoint::WilsonStable {
                        self.ip_data = &IP_HEX_WILSON_STABLE_9
                    } else {
                        self.ip_data = &IP_HEX_WILSON_CORNER_9
                    }
                }
                14 => self.ip_data = &IP_HEX_IRONS_14,
                27 => self.ip_data = &IP_HEX_LEGENDRE_27,
                _ => return Err("number of integration points is not available for Hex class"),
            },
        }
        Ok(())
    }

    /// Performs the integration of a function over the Shape's domain
    ///
    /// # General case (geo_ndim == space_ndim)
    ///
    /// ```text
    ///       ⌠   →       ⌠   → →            →
    /// res = │ f(x) dΩ = │ f(x(ξ)) ⋅ det(J)(ξ) dΩ
    ///       ⌡           ⌡
    ///       Ω           Ωref
    /// ```
    ///
    /// which is replaced by numerical integration according to:
    ///
    /// ```text
    ///       nip-1  →             →
    /// res ≈  Σ   f(ιp)) ⋅ det(J)(ιp) ⋅ wp
    ///       p=0
    /// ```
    ///
    /// where `nip` is the number of integration points, `ιp := ξp` is the reference
    /// coordinate of the integration point, and `wp` is the weight attached to the
    /// p-th integration point.
    ///
    /// # Line in multi-dimensions (geo_ndim == 1 and space_ndim > 1)
    ///
    /// ```text
    ///       ⌠               ⌠
    /// res = │ f(ell) dell = │ f(ξ(ell)) ⋅ ||Jline||(ξ) dξ
    ///       ⌡               ⌡
    ///       Ω               Ωref
    /// ```
    ///
    /// where `||Jline||` is the Euclidean norm of `Jline`.
    ///
    /// The above integral is replaced by numerical integration according to:
    ///
    /// ```text
    ///       nip-1      →               →
    /// res ≈  Σ   f(ell(ιp)) ⋅ ||Jline||(ιp) ⋅ wp
    ///       p=0
    /// ```
    ///
    /// Here, we considered a parametric coordinate `ell` which varies
    /// from 0 to `ell_max` (the length of the line) according to
    ///
    /// ```text
    ///                  ell_max
    /// ell(ξ) = (1 + ξ) ———————
    ///                     2
    ///
    ///          2 · ell
    /// ξ(ell) = ——————— - 1
    ///          ell_max
    /// ```
    ///
    /// ```text
    /// 0 ≤ ell ≤ ell_max
    ///
    /// -1 ≤ ξ ≤ +1
    /// ```
    pub fn integ(&mut self) -> Result<f64, StrError> {
        Ok(0.0)
    }

    /// Performs the integration of a function over the Shape's boundary
    ///
    /// # Boundary line in 2D (geo_ndim == 1 and space_ndim == 2)
    ///
    /// ```text
    ///       ⌠             →        ⌠           →    →
    /// res = │ q(ell) unit_n dell = │ q(ell) ⋅ (e3 × g1) dξ
    ///       ⌡                      ⌡
    ///       Γ                      Γref
    /// ```
    ///
    /// where `unit_n` is the unit normal vector.
    ///
    /// The above integral is replaced by numerical integration according to:
    ///
    /// ```text
    ///       nip-1      →       →  →     →  →
    /// res ≈  Σ   q(ell(ιp)) ⋅ (e3(ιp) × g1(ιp)) ⋅ wp
    ///       p=0
    /// ```
    ///
    /// Here, we considered a parametric coordinate `ell` which varies
    /// from 0 to `ell_max` (the length of the line) according to
    ///
    /// ```text
    ///                  ell_max
    /// ell(ξ) = (1 + ξ) ———————
    ///                     2
    ///
    ///          2 · ell
    /// ξ(ell) = ——————— - 1
    ///          ell_max
    /// ```
    ///
    /// ```text
    /// 0 ≤ ell ≤ ell_max
    ///
    /// -1 ≤ ξ ≤ +1
    /// ```
    ///
    /// # 3D surface (geo_ndim == 2 and space_ndim == 3)
    ///
    /// We can then perform integrations in the reference space as follows:
    ///
    /// ```text
    ///       ⌠   →       →      ⌠   → →      →    →
    /// res = │ f(x) unit_n dA = │ f(x(ξ)) ⋅ (g1 × g2) dξ1 dξ2
    ///       ⌡                  ⌡
    ///       Γ                  Γref
    /// ```
    ///
    /// where `unit_n` is the unit normal vector.
    ///
    /// The above integral is replaced by numerical integration according to:
    ///
    /// ```text
    ///       nip-1   →      →   →    →   →
    /// res ≈  Σ   f(ιp)) ⋅ (g1(ιp) × g2(ιp)) ⋅ wp
    ///       p=0
    /// ```
    pub fn boundary_integ(&mut self) -> Result<f64, StrError> {
        Ok(0.0)
    }
}
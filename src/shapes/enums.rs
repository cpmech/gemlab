use super::*;
use serde::{Deserialize, Serialize};

/// Defines the class of geometric shape
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Deserialize, Serialize)]
pub enum GeoClass {
    /// Lines (segments) class
    Lin,

    /// Triangles class
    Tri,

    /// Quadrilaterals class
    Qua,

    /// Tetrahedra class
    Tet,

    /// Hexahedra class
    Hex,
}

/// Defines the kind of geometric shape
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Deserialize, Serialize)]
pub enum GeoKind {
    /// Line (segment) with 2 nodes (linear functions)
    Lin2 = 1_002,

    /// Line (segment) with 3 nodes (quadratic functions)
    Lin3 = 1_003,

    /// Line (segment) with 4 nodes (cubic functions)
    Lin4 = 1_004,

    /// Line (segment) with 5 nodes (quartic functions)
    Lin5 = 1_005,

    /// Triangle with 3 nodes (linear edges)
    Tri3 = 2_003,

    /// Triangle with 6 nodes (quadratic edges)
    Tri6 = 2_006,

    /// Triangle with 10 nodes (cubic edges; interior node)
    Tri10 = 2_010,

    /// Triangle with 15 nodes (quartic edges; interior nodes)
    Tri15 = 2_015,

    /// Quadrilateral with 4 nodes (linear edges)
    Qua4 = 2_004,

    /// Quadrilateral with 8 nodes (quadratic edges)
    Qua8 = 2_008,

    /// Quadrilateral with 9 nodes (quadratic edges; interior node)
    Qua9 = 2_009,

    /// Quadrilateral with 12 nodes (cubic edges)
    Qua12 = 2_012,

    /// Quadrilateral with 16 nodes (cubic edges; interior nodes)
    Qua16 = 2_016,

    /// Quadrilateral with 17 nodes (quartic edges; interior node)
    Qua17 = 2_017,

    /// Tetrahedron with 4 nodes (linear faces)
    Tet4 = 3_004,

    /// Tetrahedron with 10 nodes (quadratic faces)
    Tet10 = 3_010,

    /// Hexahedron with 8 nodes (bilinear faces)
    Hex8 = 3_008,

    /// Hexahedron with 20 nodes (quadratic faces)
    Hex20 = 3_020,
}

impl GeoKind {
    /// Holds all enum values
    pub const VALUES: [Self; 18] = [
        // Lin
        Self::Lin2,
        Self::Lin3,
        Self::Lin4,
        Self::Lin5,
        // Tri
        Self::Tri3,
        Self::Tri6,
        Self::Tri10,
        Self::Tri15,
        // Qua
        Self::Qua4,
        Self::Qua8,
        Self::Qua9,
        Self::Qua12,
        Self::Qua16,
        Self::Qua17,
        // Tet
        Self::Tet4,
        Self::Tet10,
        // Hex
        Self::Hex8,
        Self::Hex20,
    ];
}

pub(crate) fn i32_to_fn_interp(kind: i32) -> FnInterp {
    match kind {
        // Lin
        1_002 => FnInterp(GeoKind::Lin2, Lin2::calc_interp),
        1_003 => FnInterp(GeoKind::Lin3, Lin3::calc_interp),
        1_004 => FnInterp(GeoKind::Lin4, Lin4::calc_interp),
        1_005 => FnInterp(GeoKind::Lin5, Lin5::calc_interp),
        // Tri
        2_003 => FnInterp(GeoKind::Tri3, Tri3::calc_interp),
        2_006 => FnInterp(GeoKind::Tri6, Tri6::calc_interp),
        2_010 => FnInterp(GeoKind::Tri10, Tri10::calc_interp),
        2_015 => FnInterp(GeoKind::Tri15, Tri15::calc_interp),
        // Qua
        2_004 => FnInterp(GeoKind::Qua4, Qua4::calc_interp),
        2_008 => FnInterp(GeoKind::Qua8, Qua8::calc_interp),
        2_009 => FnInterp(GeoKind::Qua9, Qua9::calc_interp),
        2_012 => FnInterp(GeoKind::Qua12, Qua12::calc_interp),
        2_016 => FnInterp(GeoKind::Qua16, Qua16::calc_interp),
        2_017 => FnInterp(GeoKind::Qua17, Qua17::calc_interp),
        // Tet
        3_004 => FnInterp(GeoKind::Tet4, Tet4::calc_interp),
        3_010 => FnInterp(GeoKind::Tet10, Tet10::calc_interp),
        // Hex
        3_008 => FnInterp(GeoKind::Hex8, Hex8::calc_interp),
        3_020 => FnInterp(GeoKind::Hex20, Hex20::calc_interp),
        _ => panic!("INTERNAL ERROR: cannot convert i32 to FnInterp"),
    }
}

pub(crate) fn i32_to_fn_deriv(kind: i32) -> FnDeriv {
    match kind {
        // Lin
        1_002 => FnDeriv(GeoKind::Lin2, Lin2::calc_deriv),
        1_003 => FnDeriv(GeoKind::Lin3, Lin3::calc_deriv),
        1_004 => FnDeriv(GeoKind::Lin4, Lin4::calc_deriv),
        1_005 => FnDeriv(GeoKind::Lin5, Lin5::calc_deriv),
        // Tri
        2_003 => FnDeriv(GeoKind::Tri3, Tri3::calc_deriv),
        2_006 => FnDeriv(GeoKind::Tri6, Tri6::calc_deriv),
        2_010 => FnDeriv(GeoKind::Tri10, Tri10::calc_deriv),
        2_015 => FnDeriv(GeoKind::Tri15, Tri15::calc_deriv),
        // Qua
        2_004 => FnDeriv(GeoKind::Qua4, Qua4::calc_deriv),
        2_008 => FnDeriv(GeoKind::Qua8, Qua8::calc_deriv),
        2_009 => FnDeriv(GeoKind::Qua9, Qua9::calc_deriv),
        2_012 => FnDeriv(GeoKind::Qua12, Qua12::calc_deriv),
        2_016 => FnDeriv(GeoKind::Qua16, Qua16::calc_deriv),
        2_017 => FnDeriv(GeoKind::Qua17, Qua17::calc_deriv),
        // Tet
        3_004 => FnDeriv(GeoKind::Tet4, Tet4::calc_deriv),
        3_010 => FnDeriv(GeoKind::Tet10, Tet10::calc_deriv),
        // Hex
        3_008 => FnDeriv(GeoKind::Hex8, Hex8::calc_deriv),
        3_020 => FnDeriv(GeoKind::Hex20, Hex20::calc_deriv),
        _ => panic!("INTERNAL ERROR: cannot convert i32 to FnDeriv"),
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{i32_to_fn_deriv, i32_to_fn_interp};
    use crate::shapes::{GeoClass, GeoKind};
    use std::collections::HashSet;

    #[test]
    fn derive_works() {
        let class = GeoClass::Tri.clone();
        let kind = GeoKind::Tri6.clone();
        let classes = HashSet::from([GeoClass::Tri, GeoClass::Qua]);
        let kinds = HashSet::from([GeoKind::Tri3, GeoKind::Qua4]);
        assert_eq!(class, GeoClass::Tri);
        assert_eq!(kind, GeoKind::Tri6);
        assert_eq!(format!("{:?}", class), "Tri");
        assert_eq!(format!("{:?}", kind), "Tri6");
        assert_eq!(classes.contains(&GeoClass::Tri), true);
        assert_eq!(kinds.contains(&GeoKind::Tri3), true);
    }

    #[test]
    fn data_is_consistent() {
        for kind in GeoKind::VALUES {
            match kind {
                // Lin
                GeoKind::Lin2 => {
                    let k = kind as i32;
                    let nnode = k - 1000;
                    assert_eq!(nnode, 2);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Lin2);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Lin2);
                }
                GeoKind::Lin3 => {
                    let k = kind as i32;
                    let nnode = k - 1000;
                    assert_eq!(nnode, 3);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Lin3);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Lin3);
                }
                GeoKind::Lin4 => {
                    let k = kind as i32;
                    let nnode = k - 1000;
                    assert_eq!(nnode, 4);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Lin4);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Lin4);
                }
                GeoKind::Lin5 => {
                    let k = kind as i32;
                    let nnode = k - 1000;
                    assert_eq!(nnode, 5);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Lin5);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Lin5);
                }

                // Tri
                GeoKind::Tri3 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 3);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tri3);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tri3);
                }
                GeoKind::Tri6 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 6);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tri6);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tri6);
                }
                GeoKind::Tri10 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 10);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tri10);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tri10);
                }
                GeoKind::Tri15 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 15);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tri15);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tri15);
                }

                // Qua
                GeoKind::Qua4 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 4);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua4);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua4);
                }
                GeoKind::Qua8 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 8);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua8);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua8);
                }
                GeoKind::Qua9 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 9);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua9);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua9);
                }
                GeoKind::Qua12 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 12);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua12);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua12);
                }
                GeoKind::Qua16 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 16);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua16);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua16);
                }
                GeoKind::Qua17 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 17);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua17);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua17);
                }

                // Tet
                GeoKind::Tet4 => {
                    let k = kind as i32;
                    let nnode = k - 3000;
                    assert_eq!(nnode, 4);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tet4);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tet4);
                }
                GeoKind::Tet10 => {
                    let k = kind as i32;
                    let nnode = k - 3000;
                    assert_eq!(nnode, 10);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tet10);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tet10);
                }

                // Hex
                GeoKind::Hex8 => {
                    let k = kind as i32;
                    let nnode = k - 3000;
                    assert_eq!(nnode, 8);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Hex8);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Hex8);
                }
                GeoKind::Hex20 => {
                    let k = kind as i32;
                    let nnode = k - 3000;
                    assert_eq!(nnode, 20);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Hex20);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Hex20);
                }
            }
        }
    }
}

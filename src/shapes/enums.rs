use super::*;
use crate::StrError;
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

/// Returns the GeoClass and GeoKind for given geo_ndim and nnode
pub fn geo_class_and_kind(geo_ndim: usize, nnode: usize) -> Result<(GeoClass, GeoKind), StrError> {
    match (geo_ndim, nnode) {
        // Lin
        (1, 2) => Ok((GeoClass::Lin, GeoKind::Lin2)),
        (1, 3) => Ok((GeoClass::Lin, GeoKind::Lin3)),
        (1, 4) => Ok((GeoClass::Lin, GeoKind::Lin4)),
        (1, 5) => Ok((GeoClass::Lin, GeoKind::Lin5)),

        // Tri
        (2, 3) => Ok((GeoClass::Tri, GeoKind::Tri3)),
        (2, 6) => Ok((GeoClass::Tri, GeoKind::Tri6)),
        (2, 10) => Ok((GeoClass::Tri, GeoKind::Tri10)),
        (2, 15) => Ok((GeoClass::Tri, GeoKind::Tri15)),

        // Qua
        (2, 4) => Ok((GeoClass::Qua, GeoKind::Qua4)),
        (2, 8) => Ok((GeoClass::Qua, GeoKind::Qua8)),
        (2, 9) => Ok((GeoClass::Qua, GeoKind::Qua9)),
        (2, 12) => Ok((GeoClass::Qua, GeoKind::Qua12)),
        (2, 16) => Ok((GeoClass::Qua, GeoKind::Qua16)),
        (2, 17) => Ok((GeoClass::Qua, GeoKind::Qua17)),

        // Tet
        (3, 4) => Ok((GeoClass::Tet, GeoKind::Tet4)),
        (3, 10) => Ok((GeoClass::Tet, GeoKind::Tet10)),

        // Hex
        (3, 8) => Ok((GeoClass::Hex, GeoKind::Hex8)),
        (3, 20) => Ok((GeoClass::Hex, GeoKind::Hex20)),
        _ => Err("(geo_ndim,nnode) combination is invalid"),
    }
}

/// Converts i32 to FnInterp
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

/// Converts i32 to FnDeriv
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
    use crate::{
        shapes::{geo_class_and_kind, GeoClass, GeoKind},
        StrError,
    };
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
    #[should_panic(expected = "INTERNAL ERROR: cannot convert i32 to FnInterp")]
    fn i32_to_fn_interp_panics_on_wrong_input() {
        i32_to_fn_interp(-1);
    }

    #[test]
    #[should_panic(expected = "INTERNAL ERROR: cannot convert i32 to FnDeriv")]
    fn i32_to_fn_deriv_panics_on_wrong_input() {
        i32_to_fn_deriv(-1);
    }

    #[test]
    fn data_is_consistent() -> Result<(), StrError> {
        for kind in GeoKind::VALUES {
            match kind {
                // Lin
                GeoKind::Lin2 => {
                    let k = kind as i32;
                    let nnode = k - 1000;
                    assert_eq!(nnode, 2);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Lin2);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Lin2);
                    assert_eq!(geo_class_and_kind(1, 2)?, (GeoClass::Lin, kind));
                }
                GeoKind::Lin3 => {
                    let k = kind as i32;
                    let nnode = k - 1000;
                    assert_eq!(nnode, 3);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Lin3);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Lin3);
                    assert_eq!(geo_class_and_kind(1, 3)?, (GeoClass::Lin, kind));
                }
                GeoKind::Lin4 => {
                    let k = kind as i32;
                    let nnode = k - 1000;
                    assert_eq!(nnode, 4);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Lin4);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Lin4);
                    assert_eq!(geo_class_and_kind(1, 4)?, (GeoClass::Lin, kind));
                }
                GeoKind::Lin5 => {
                    let k = kind as i32;
                    let nnode = k - 1000;
                    assert_eq!(nnode, 5);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Lin5);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Lin5);
                    assert_eq!(geo_class_and_kind(1, 5)?, (GeoClass::Lin, kind));
                }

                // Tri
                GeoKind::Tri3 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 3);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tri3);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tri3);
                    assert_eq!(geo_class_and_kind(2, 3)?, (GeoClass::Tri, kind));
                }
                GeoKind::Tri6 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 6);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tri6);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tri6);
                    assert_eq!(geo_class_and_kind(2, 6)?, (GeoClass::Tri, kind));
                }
                GeoKind::Tri10 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 10);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tri10);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tri10);
                    assert_eq!(geo_class_and_kind(2, 10)?, (GeoClass::Tri, kind));
                }
                GeoKind::Tri15 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 15);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tri15);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tri15);
                    assert_eq!(geo_class_and_kind(2, 15)?, (GeoClass::Tri, kind));
                }

                // Qua
                GeoKind::Qua4 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 4);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua4);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua4);
                    assert_eq!(geo_class_and_kind(2, 4)?, (GeoClass::Qua, kind));
                }
                GeoKind::Qua8 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 8);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua8);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua8);
                    assert_eq!(geo_class_and_kind(2, 8)?, (GeoClass::Qua, kind));
                }
                GeoKind::Qua9 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 9);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua9);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua9);
                    assert_eq!(geo_class_and_kind(2, 9)?, (GeoClass::Qua, kind));
                }
                GeoKind::Qua12 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 12);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua12);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua12);
                    assert_eq!(geo_class_and_kind(2, 12)?, (GeoClass::Qua, kind));
                }
                GeoKind::Qua16 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 16);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua16);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua16);
                    assert_eq!(geo_class_and_kind(2, 16)?, (GeoClass::Qua, kind));
                }
                GeoKind::Qua17 => {
                    let k = kind as i32;
                    let nnode = k - 2000;
                    assert_eq!(nnode, 17);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Qua17);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Qua17);
                    assert_eq!(geo_class_and_kind(2, 17)?, (GeoClass::Qua, kind));
                }

                // Tet
                GeoKind::Tet4 => {
                    let k = kind as i32;
                    let nnode = k - 3000;
                    assert_eq!(nnode, 4);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tet4);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tet4);
                    assert_eq!(geo_class_and_kind(3, 4)?, (GeoClass::Tet, kind));
                }
                GeoKind::Tet10 => {
                    let k = kind as i32;
                    let nnode = k - 3000;
                    assert_eq!(nnode, 10);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Tet10);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Tet10);
                    assert_eq!(geo_class_and_kind(3, 10)?, (GeoClass::Tet, kind));
                }

                // Hex
                GeoKind::Hex8 => {
                    let k = kind as i32;
                    let nnode = k - 3000;
                    assert_eq!(nnode, 8);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Hex8);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Hex8);
                    assert_eq!(geo_class_and_kind(3, 8)?, (GeoClass::Hex, kind));
                }
                GeoKind::Hex20 => {
                    let k = kind as i32;
                    let nnode = k - 3000;
                    assert_eq!(nnode, 20);
                    assert_eq!(i32_to_fn_interp(k).0, GeoKind::Hex20);
                    assert_eq!(i32_to_fn_deriv(k).0, GeoKind::Hex20);
                    assert_eq!(geo_class_and_kind(3, 20)?, (GeoClass::Hex, kind));
                }
            }
        }
        Ok(())
    }
}

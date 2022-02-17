/// Defines the class of geometric shape
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
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
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum GeoKind {
    /// Line (segment) with 2 nodes (linear functions)
    Lin2,

    /// Line (segment) with 3 nodes (quadratic functions)
    Lin3,

    /// Line (segment) with 4 nodes (cubic functions)
    Lin4,

    /// Line (segment) with 5 nodes (quartic functions)
    Lin5,

    /// Triangle with 3 nodes (linear edges)
    Tri3,

    /// Triangle with 6 nodes (quadratic edges)
    Tri6,

    /// Triangle with 10 nodes (cubic edges; interior node)
    Tri10,

    /// Triangle with 15 nodes (quartic edges; interior nodes)
    Tri15,

    /// Quadrilateral with 4 nodes (linear edges)
    Qua4,

    /// Quadrilateral with 8 nodes (quadratic edges)
    Qua8,

    /// Quadrilateral with 9 nodes (quadratic edges; interior node)
    Qua9,

    /// Quadrilateral with 12 nodes (cubic edges)
    Qua12,

    /// Quadrilateral with 16 nodes (cubic edges; interior nodes)
    Qua16,

    /// Quadrilateral with 17 nodes (quartic edges; interior node)
    Qua17,

    /// Tetrahedron with 4 nodes (linear faces)
    Tet4,

    /// Tetrahedron with 10 nodes (quadratic faces)
    Tet10,

    /// Hexahedron with 8 nodes (bilinear faces)
    Hex8,

    /// Hexahedron with 20 nodes (quadratic faces)
    Hex20,
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
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
}

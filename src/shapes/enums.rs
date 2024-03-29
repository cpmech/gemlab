use super::{
    Hex20, Hex32, Hex8, Lin2, Lin3, Lin4, Lin5, Qua12, Qua16, Qua17, Qua4, Qua8, Qua9, Tet10, Tet20, Tet4, Tri10,
    Tri15, Tri3, Tri6,
};
use crate::StrError;
use russell_lab::{Matrix, Vector};
use serde::{Deserialize, Serialize};

/// Defines a function to calculate interpolation functions
///
/// # Input
///
/// * `interp` -- vector with the resulting interpolation functions; its length mut be equal to `nnode`
/// * `ksi` -- set of coordinates in the reference space; its length must be ≥ `geo_ndim`
///
/// # Panics
///
/// None of the dimensions are checked, thus a panic my occur if the
/// arguments have dimensions incompatible with the related [GeoKind].
pub type FnInterp = fn(interp: &mut Vector, ksi: &[f64]);

/// Defines a function to calculate the derivatives of interpolation functions
///
/// # Input
///
/// * `deriv` -- matrix with the resulting set of derivatives of interpolation functions
///   with respect to the reference (natural) coordinates; its dimensions are `(nnode,geo_ndim)`
/// * `ksi` -- set of coordinates in the reference space; its length must be ≥ `geo_ndim`
///
/// # Panics
///
/// None of the dimensions are checked, thus a panic my occur if the
/// arguments have dimensions incompatible with the related [GeoKind].
pub type FnDeriv = fn(deriv: &mut Matrix, ksi: &[f64]);

/// Defines the geometry case regarding the number of dimensions (geo vs space) (Cable, Shell, Solid)
///
/// The following table shows what combinations of geometry-number-of-dimensions (`geo_ndim`) and
/// space-number-of-dimensions (`space_ndim`) are possible. There are three cases:
///
/// 1. Case `CABLE` -- `geo_ndim = 1` and `space_ndim = 2 or 3`; e.g., line in 2D or 3D (cables and rods)
/// 2. Case `SHELL` -- `geo_ndim = 2` and `space_ndim = 3`; e.g. Tri or Qua in 3D (shells and surfaces)
/// 3. Case `SOLID` -- `geo_ndim = space_ndim`; e.g., Tri and Qua in 2D or Tet and Hex in 3D
///
/// | `geo_ndim` | `space_ndim = 2` | `space_ndim = 3` |
/// |:----------:|:----------------:|:----------------:|
/// |     1      |     `CABLE`      |     `CABLE`      |
/// |     2      |     `SOLID`      |     `SHELL`      |
/// |     3      |    impossible    |     `SOLID`      |
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GeoCase {
    /// Represents cables (rods) in 2D or 3D with geo_ndim = 1 and space_ndim = 2 or 3
    Cable,

    /// Represents shells (surfaces) in 3D with geo_ndim = 2 and space_ndim = 3
    Shell,

    /// Represents solids in 2D or 3D with geo_ndim = space_ndim
    Solid,
}

/// Returns the geometry case given the geo and space dimensions
///
/// # Panics
///
/// 1. `space_ndim` must be 2 or 3; otherwise a panic will occur
/// 2. This function will panic if `geo_ndim > space_ndim` (impossible case)
#[inline]
pub fn geo_case(geo_ndim: usize, space_ndim: usize) -> GeoCase {
    assert!(space_ndim >= 2 && space_ndim <= 3);
    assert!(geo_ndim <= space_ndim);
    if geo_ndim == space_ndim {
        GeoCase::Solid
    } else if geo_ndim == 1 {
        GeoCase::Cable
    } else {
        GeoCase::Shell
    }
}

/// Defines the class of the geometric shape (Lin, Tri, Qua, Tet, Hex)
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

/// Defines the kind of the geometric shape (Lin2, ... Tri3, ..., Qua4, ... Tet4, ..., Hex8, ...)
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Deserialize, Serialize)]
pub enum GeoKind {
    /// Line (segment) with 2 nodes (linear functions)
    Lin2,

    /// Line (segment) with 3 nodes (quadratic functions; with 1 interior/middle node)
    Lin3,

    /// Line (segment) with 4 nodes (cubic functions; with 2 interior nodes)
    Lin4,

    /// Line (segment) with 5 nodes (quartic functions; with 3 interior nodes)
    Lin5,

    /// Triangle with 3 nodes (linear edges)
    Tri3,

    /// Triangle with 6 nodes (quadratic edges)
    Tri6,

    /// Triangle with 10 nodes (cubic edges; with 1 interior node)
    Tri10,

    /// Triangle with 15 nodes (quartic edges; with 3 interior nodes)
    Tri15,

    /// Quadrilateral with 4 nodes (linear edges)
    Qua4,

    /// Quadrilateral with 8 nodes (quadratic edges)
    Qua8,

    /// Quadrilateral with 9 nodes (quadratic edges; with 1 interior node)
    Qua9,

    /// Quadrilateral with 12 nodes (cubic edges)
    Qua12,

    /// Quadrilateral with 16 nodes (cubic edges; with 4 interior nodes)
    Qua16,

    /// Quadrilateral with 17 nodes (quartic edges; with 1 interior node)
    Qua17,

    /// Tetrahedron with 4 nodes (linear faces = Tri3)
    Tet4,

    /// Tetrahedron with 10 nodes (quadratic faces = Tri6)
    Tet10,

    /// Tetrahedron with 20 nodes (cubic faces = Tri10)
    Tet20,

    /// Hexahedron with 8 nodes (bilinear faces = Qua4)
    Hex8,

    /// Hexahedron with 20 nodes (quadratic faces = Qua8)
    Hex20,

    /// Hexahedron with 32 nodes (cubic faces = Qua12)
    Hex32,
}

impl GeoKind {
    /// Returns a new GeoKind given a string representation
    ///
    /// **Important:** The input string must be **lowercase**.
    ///
    /// **Lin**
    ///
    /// * `lin2` -- Line with 2 nodes (linear functions)
    /// * `lin3` -- Line with 2 nodes (quadratic functions)
    /// * `lin4` -- Line with 2 nodes (cubic functions)
    /// * `lin5` -- Line with 2 nodes (quartic functions)
    ///
    /// **Tri**
    ///
    /// * `tri3` -- Triangle with 3 nodes (linear edges)
    /// * `tri6` -- Triangle with 6 nodes (quadratic edges)
    /// * `tri10` -- Triangle with 10 nodes (cubic edges; with 1 interior node)
    /// * `tri15` -- Triangle with 15 nodes (quartic edges; with 3 interior nodes)
    ///
    /// **Qua**
    ///
    /// * `qua4` -- Quadrilateral with 4 nodes (linear edges)
    /// * `qua8` -- Quadrilateral with 8 nodes (quadratic edges)
    /// * `qua9` -- Quadrilateral with 9 nodes (quadratic edges; with 1 interior node)
    /// * `qua12` -- Quadrilateral with 12 nodes (cubic edges)
    /// * `qua16` -- Quadrilateral with 16 nodes (cubic edges; with 4 interior nodes)
    /// * `qua17` -- Quadrilateral with 17 nodes (quartic edges; with 1 interior node)
    ///
    /// **Tet**
    ///
    /// * `tet4` -- Tetrahedron with 4 nodes (bilinear faces = tri3)
    /// * `tet10` -- Tetrahedron with 10 nodes (quadratic faces = tri6)
    /// * `tet20` -- Tetrahedron with 20 nodes (cubic faces = tri10)
    ///
    /// **Hex**
    ///
    /// * `hex8` -- Hexahedron with 8 nodes (bilinear faces = qua4)
    /// * `hex20` -- Hexahedron with 20 nodes (quadratic faces = qua8)
    /// * `hex32` -- Hexahedron with 32 nodes (cubic faces = qua12)
    pub fn from(kind: &str) -> Result<Self, StrError> {
        match kind {
            // Lin
            "lin2" => Ok(Self::Lin2),
            "lin3" => Ok(Self::Lin3),
            "lin4" => Ok(Self::Lin4),
            "lin5" => Ok(Self::Lin5),
            // Tri
            "tri3" => Ok(Self::Tri3),
            "tri6" => Ok(Self::Tri6),
            "tri10" => Ok(Self::Tri10),
            "tri15" => Ok(Self::Tri15),
            // Qua
            "qua4" => Ok(Self::Qua4),
            "qua8" => Ok(Self::Qua8),
            "qua9" => Ok(Self::Qua9),
            "qua12" => Ok(Self::Qua12),
            "qua16" => Ok(Self::Qua16),
            "qua17" => Ok(Self::Qua17),
            // Tet
            "tet4" => Ok(Self::Tet4),
            "tet10" => Ok(Self::Tet10),
            "tet20" => Ok(Self::Tet20),
            // Hex
            "hex8" => Ok(Self::Hex8),
            "hex20" => Ok(Self::Hex20),
            "hex32" => Ok(Self::Hex32),
            _ => Err("string representation of GeoKind is incorrect"),
        }
    }

    /// Returns the string representation
    pub fn to_string(&self) -> String {
        match self {
            // Lin
            Self::Lin2 => "lin2".to_string(),
            Self::Lin3 => "lin3".to_string(),
            Self::Lin4 => "lin4".to_string(),
            Self::Lin5 => "lin5".to_string(),
            // Tri
            Self::Tri3 => "tri3".to_string(),
            Self::Tri6 => "tri6".to_string(),
            Self::Tri10 => "tri10".to_string(),
            Self::Tri15 => "tri15".to_string(),
            // Qua
            Self::Qua4 => "qua4".to_string(),
            Self::Qua8 => "qua8".to_string(),
            Self::Qua9 => "qua9".to_string(),
            Self::Qua12 => "qua12".to_string(),
            Self::Qua16 => "qua16".to_string(),
            Self::Qua17 => "qua17".to_string(),
            // Tet
            Self::Tet4 => "tet4".to_string(),
            Self::Tet10 => "tet10".to_string(),
            Self::Tet20 => "tet20".to_string(),
            // Hex
            Self::Hex8 => "hex8".to_string(),
            Self::Hex20 => "hex20".to_string(),
            Self::Hex32 => "hex32".to_string(),
        }
    }

    /// Returns the class
    pub fn class(&self) -> GeoClass {
        match self {
            // Lin
            Self::Lin2 => GeoClass::Lin,
            Self::Lin3 => GeoClass::Lin,
            Self::Lin4 => GeoClass::Lin,
            Self::Lin5 => GeoClass::Lin,
            // Tri
            Self::Tri3 => GeoClass::Tri,
            Self::Tri6 => GeoClass::Tri,
            Self::Tri10 => GeoClass::Tri,
            Self::Tri15 => GeoClass::Tri,
            // Qua
            Self::Qua4 => GeoClass::Qua,
            Self::Qua8 => GeoClass::Qua,
            Self::Qua9 => GeoClass::Qua,
            Self::Qua12 => GeoClass::Qua,
            Self::Qua16 => GeoClass::Qua,
            Self::Qua17 => GeoClass::Qua,
            // Tet
            Self::Tet4 => GeoClass::Tet,
            Self::Tet10 => GeoClass::Tet,
            Self::Tet20 => GeoClass::Tet,
            // Hex
            Self::Hex8 => GeoClass::Hex,
            Self::Hex20 => GeoClass::Hex,
            Self::Hex32 => GeoClass::Hex,
        }
    }

    /// Returns the (geometry) ndim
    pub fn ndim(&self) -> usize {
        match self {
            // Lin
            Self::Lin2 => Lin2::GEO_NDIM,
            Self::Lin3 => Lin3::GEO_NDIM,
            Self::Lin4 => Lin4::GEO_NDIM,
            Self::Lin5 => Lin5::GEO_NDIM,
            // Tri
            Self::Tri3 => Tri3::GEO_NDIM,
            Self::Tri6 => Tri6::GEO_NDIM,
            Self::Tri10 => Tri10::GEO_NDIM,
            Self::Tri15 => Tri15::GEO_NDIM,
            // Qua
            Self::Qua4 => Qua4::GEO_NDIM,
            Self::Qua8 => Qua8::GEO_NDIM,
            Self::Qua9 => Qua9::GEO_NDIM,
            Self::Qua12 => Qua12::GEO_NDIM,
            Self::Qua16 => Qua16::GEO_NDIM,
            Self::Qua17 => Qua17::GEO_NDIM,
            // Tet
            Self::Tet4 => Tet4::GEO_NDIM,
            Self::Tet10 => Tet10::GEO_NDIM,
            Self::Tet20 => Tet20::GEO_NDIM,
            // Hex
            Self::Hex8 => Hex8::GEO_NDIM,
            Self::Hex20 => Hex20::GEO_NDIM,
            Self::Hex32 => Hex32::GEO_NDIM,
        }
    }

    /// Returns the number of nodes
    pub fn nnode(&self) -> usize {
        match self {
            // Lin
            Self::Lin2 => Lin2::NNODE,
            Self::Lin3 => Lin3::NNODE,
            Self::Lin4 => Lin4::NNODE,
            Self::Lin5 => Lin5::NNODE,
            // Tri
            Self::Tri3 => Tri3::NNODE,
            Self::Tri6 => Tri6::NNODE,
            Self::Tri10 => Tri10::NNODE,
            Self::Tri15 => Tri15::NNODE,
            // Qua
            Self::Qua4 => Qua4::NNODE,
            Self::Qua8 => Qua8::NNODE,
            Self::Qua9 => Qua9::NNODE,
            Self::Qua12 => Qua12::NNODE,
            Self::Qua16 => Qua16::NNODE,
            Self::Qua17 => Qua17::NNODE,
            // Tet
            Self::Tet4 => Tet4::NNODE,
            Self::Tet10 => Tet10::NNODE,
            Self::Tet20 => Tet20::NNODE,
            // Hex
            Self::Hex8 => Hex8::NNODE,
            Self::Hex20 => Hex20::NNODE,
            Self::Hex32 => Hex32::NNODE,
        }
    }

    /// Returns the number of edges
    pub fn nedge(&self) -> usize {
        match self {
            // Lin
            Self::Lin2 => Lin2::NEDGE,
            Self::Lin3 => Lin3::NEDGE,
            Self::Lin4 => Lin4::NEDGE,
            Self::Lin5 => Lin5::NEDGE,
            // Tri
            Self::Tri3 => Tri3::NEDGE,
            Self::Tri6 => Tri6::NEDGE,
            Self::Tri10 => Tri10::NEDGE,
            Self::Tri15 => Tri15::NEDGE,
            // Qua
            Self::Qua4 => Qua4::NEDGE,
            Self::Qua8 => Qua8::NEDGE,
            Self::Qua9 => Qua9::NEDGE,
            Self::Qua12 => Qua12::NEDGE,
            Self::Qua16 => Qua16::NEDGE,
            Self::Qua17 => Qua17::NEDGE,
            // Tet
            Self::Tet4 => Tet4::NEDGE,
            Self::Tet10 => Tet10::NEDGE,
            Self::Tet20 => Tet20::NEDGE,
            // Hex
            Self::Hex8 => Hex8::NEDGE,
            Self::Hex20 => Hex20::NEDGE,
            Self::Hex32 => Hex32::NEDGE,
        }
    }

    /// Returns the number of faces
    pub fn nface(&self) -> usize {
        match self {
            // Lin
            Self::Lin2 => Lin2::NFACE,
            Self::Lin3 => Lin3::NFACE,
            Self::Lin4 => Lin4::NFACE,
            Self::Lin5 => Lin5::NFACE,
            // Tri
            Self::Tri3 => Tri3::NFACE,
            Self::Tri6 => Tri6::NFACE,
            Self::Tri10 => Tri10::NFACE,
            Self::Tri15 => Tri15::NFACE,
            // Qua
            Self::Qua4 => Qua4::NFACE,
            Self::Qua8 => Qua8::NFACE,
            Self::Qua9 => Qua9::NFACE,
            Self::Qua12 => Qua12::NFACE,
            Self::Qua16 => Qua16::NFACE,
            Self::Qua17 => Qua17::NFACE,
            // Tet
            Self::Tet4 => Tet4::NFACE,
            Self::Tet10 => Tet10::NFACE,
            Self::Tet20 => Tet20::NFACE,
            // Hex
            Self::Hex8 => Hex8::NFACE,
            Self::Hex20 => Hex20::NFACE,
            Self::Hex32 => Hex32::NFACE,
        }
    }

    /// Returns the number of nodes on the edge
    pub fn edge_nnode(&self) -> usize {
        match self {
            // Lin
            Self::Lin2 => Lin2::EDGE_NNODE,
            Self::Lin3 => Lin3::EDGE_NNODE,
            Self::Lin4 => Lin4::EDGE_NNODE,
            Self::Lin5 => Lin5::EDGE_NNODE,
            // Tri
            Self::Tri3 => Tri3::EDGE_NNODE,
            Self::Tri6 => Tri6::EDGE_NNODE,
            Self::Tri10 => Tri10::EDGE_NNODE,
            Self::Tri15 => Tri15::EDGE_NNODE,
            // Qua
            Self::Qua4 => Qua4::EDGE_NNODE,
            Self::Qua8 => Qua8::EDGE_NNODE,
            Self::Qua9 => Qua9::EDGE_NNODE,
            Self::Qua12 => Qua12::EDGE_NNODE,
            Self::Qua16 => Qua16::EDGE_NNODE,
            Self::Qua17 => Qua17::EDGE_NNODE,
            // Tet
            Self::Tet4 => Tet4::EDGE_NNODE,
            Self::Tet10 => Tet10::EDGE_NNODE,
            Self::Tet20 => Tet20::EDGE_NNODE,
            // Hex
            Self::Hex8 => Hex8::EDGE_NNODE,
            Self::Hex20 => Hex20::EDGE_NNODE,
            Self::Hex32 => Hex32::EDGE_NNODE,
        }
    }

    /// Returns the number of nodes on the face
    pub fn face_nnode(&self) -> usize {
        match self {
            // Lin
            Self::Lin2 => Lin2::FACE_NNODE,
            Self::Lin3 => Lin3::FACE_NNODE,
            Self::Lin4 => Lin4::FACE_NNODE,
            Self::Lin5 => Lin5::FACE_NNODE,
            // Tri
            Self::Tri3 => Tri3::FACE_NNODE,
            Self::Tri6 => Tri6::FACE_NNODE,
            Self::Tri10 => Tri10::FACE_NNODE,
            Self::Tri15 => Tri15::FACE_NNODE,
            // Qua
            Self::Qua4 => Qua4::FACE_NNODE,
            Self::Qua8 => Qua8::FACE_NNODE,
            Self::Qua9 => Qua9::FACE_NNODE,
            Self::Qua12 => Qua12::FACE_NNODE,
            Self::Qua16 => Qua16::FACE_NNODE,
            Self::Qua17 => Qua17::FACE_NNODE,
            // Tet
            Self::Tet4 => Tet4::FACE_NNODE,
            Self::Tet10 => Tet10::FACE_NNODE,
            Self::Tet20 => Tet20::FACE_NNODE,
            // Hex
            Self::Hex8 => Hex8::FACE_NNODE,
            Self::Hex20 => Hex20::FACE_NNODE,
            Self::Hex32 => Hex32::FACE_NNODE,
        }
    }

    /// Returns the number of edges on the face
    pub fn face_nedge(&self) -> usize {
        match self {
            // Lin
            Self::Lin2 => Lin2::FACE_NEDGE,
            Self::Lin3 => Lin3::FACE_NEDGE,
            Self::Lin4 => Lin4::FACE_NEDGE,
            Self::Lin5 => Lin5::FACE_NEDGE,
            // Tri
            Self::Tri3 => Tri3::FACE_NEDGE,
            Self::Tri6 => Tri6::FACE_NEDGE,
            Self::Tri10 => Tri10::FACE_NEDGE,
            Self::Tri15 => Tri15::FACE_NEDGE,
            // Qua
            Self::Qua4 => Qua4::FACE_NEDGE,
            Self::Qua8 => Qua8::FACE_NEDGE,
            Self::Qua9 => Qua9::FACE_NEDGE,
            Self::Qua12 => Qua12::FACE_NEDGE,
            Self::Qua16 => Qua16::FACE_NEDGE,
            Self::Qua17 => Qua17::FACE_NEDGE,
            // Tet
            Self::Tet4 => Tet4::FACE_NEDGE,
            Self::Tet10 => Tet10::FACE_NEDGE,
            Self::Tet20 => Tet20::FACE_NEDGE,
            // Hex
            Self::Hex8 => Hex8::FACE_NEDGE,
            Self::Hex20 => Hex20::FACE_NEDGE,
            Self::Hex32 => Hex32::FACE_NEDGE,
        }
    }

    /// Returns the kind of edges
    pub fn edge_kind(&self) -> Option<GeoKind> {
        match self {
            // Lin
            Self::Lin2 => None,
            Self::Lin3 => None,
            Self::Lin4 => None,
            Self::Lin5 => None,
            // Tri
            Self::Tri3 => Some(GeoKind::Lin2),
            Self::Tri6 => Some(GeoKind::Lin3),
            Self::Tri10 => Some(GeoKind::Lin4),
            Self::Tri15 => Some(GeoKind::Lin5),
            // Qua
            Self::Qua4 => Some(GeoKind::Lin2),
            Self::Qua8 => Some(GeoKind::Lin3),
            Self::Qua9 => Some(GeoKind::Lin3),
            Self::Qua12 => Some(GeoKind::Lin4),
            Self::Qua16 => Some(GeoKind::Lin4),
            Self::Qua17 => Some(GeoKind::Lin5),
            // Tet
            Self::Tet4 => Some(GeoKind::Lin2),
            Self::Tet10 => Some(GeoKind::Lin3),
            Self::Tet20 => Some(GeoKind::Lin4),
            // Hex
            Self::Hex8 => Some(GeoKind::Lin2),
            Self::Hex20 => Some(GeoKind::Lin3),
            Self::Hex32 => Some(GeoKind::Lin4),
        }
    }

    /// Returns the kind of faces
    pub fn face_kind(&self) -> Option<GeoKind> {
        match self {
            // Lin
            Self::Lin2 => None,
            Self::Lin3 => None,
            Self::Lin4 => None,
            Self::Lin5 => None,
            // Tri
            Self::Tri3 => None,
            Self::Tri6 => None,
            Self::Tri10 => None,
            Self::Tri15 => None,
            // Qua
            Self::Qua4 => None,
            Self::Qua8 => None,
            Self::Qua9 => None,
            Self::Qua12 => None,
            Self::Qua16 => None,
            Self::Qua17 => None,
            // Tet
            Self::Tet4 => Some(GeoKind::Tri3),
            Self::Tet10 => Some(GeoKind::Tri6),
            Self::Tet20 => Some(GeoKind::Tri10),
            // Hex
            Self::Hex8 => Some(GeoKind::Qua4),
            Self::Hex20 => Some(GeoKind::Qua8),
            Self::Hex32 => Some(GeoKind::Qua12),
        }
    }

    /// Returns the compatible lower-order version (e.g., Tri6 => Tri3, Qua17 => Qua9)
    pub fn lower_order(&self) -> Option<GeoKind> {
        match self {
            // Lin
            Self::Lin2 => None,
            Self::Lin3 => Some(GeoKind::Lin2),
            Self::Lin4 => None,
            Self::Lin5 => Some(GeoKind::Lin3),
            // Tri
            Self::Tri3 => None,
            Self::Tri6 => Some(GeoKind::Tri3),
            Self::Tri10 => None,
            Self::Tri15 => Some(GeoKind::Tri6),
            // Qua
            Self::Qua4 => None,
            Self::Qua8 => Some(GeoKind::Qua4),
            Self::Qua9 => Some(GeoKind::Qua4),
            Self::Qua12 => None,
            Self::Qua16 => None,
            Self::Qua17 => Some(GeoKind::Qua9),
            // Tet
            Self::Tet4 => None,
            Self::Tet10 => Some(GeoKind::Tet4),
            Self::Tet20 => None,
            // Hex
            Self::Hex8 => None,
            Self::Hex20 => Some(GeoKind::Hex8),
            Self::Hex32 => None,
        }
    }

    /// Returns the local id of node on edge (outward normals)
    ///
    /// # Input
    ///
    /// * `e` -- index of edge in [0, nedge-1]
    /// * `i` -- index of local node [0, edge_nnode-1]
    ///
    /// **Note:** For 2D GeoKinds, the sequence of edge nodes is such that the normals are **outward**.
    /// For 3D GeoKinds, the sequence just follows the associate Lin GeoKind.
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::GeoKind;
    ///
    /// //         2            OUTWARD NORMALS
    /// //   3-----6-----2          p0 p1 p2
    /// //   |           |      e:0 [1, 0, 4]
    /// //   |           |      e:1 [2, 1, 5]
    /// // 3 7           5 1    e:2 [3, 2, 6]
    /// //   |           |      e:3 [0, 3, 7]
    /// //   |           |
    /// //   0-----4-----1
    /// //         0
    /// let kind = GeoKind::Qua8;
    /// let edges: Vec<_> = (0..kind.nedge())
    ///     .map(|e| (kind.edge_node_id(e, 0),
    ///               kind.edge_node_id(e, 1),
    ///               kind.edge_node_id(e, 2)))
    ///     .collect();
    /// assert_eq!(
    ///     edges,
    ///     &[(1, 0, 4), (2, 1, 5), (3, 2, 6), (0, 3, 7)]
    /// );
    /// ```
    ///
    /// ```
    /// use gemlab::shapes::GeoKind;
    ///
    /// //          .4-----[7]------7
    /// //        ,' |            ,'|
    /// //      [4] [8]         [6] |  [#] indicates local
    /// //    ,'     |        ,'    |      edge number "e"
    /// //  5'========[5]===6'     [11]
    /// //  |               |       |
    /// //  |        |      |       |
    /// // [9]      ,0-[3]- | - - - 3
    /// //  |     ,'       [10]  [2]
    /// //  |   [0]         |   ,'
    /// //  | ,'            | ,'
    /// //  1'-----[1]------2'
    /// let kind = GeoKind::Hex8;
    /// let edges: Vec<_> = (0..kind.nedge())
    ///     .map(|e| (kind.edge_node_id(e, 0), kind.edge_node_id(e, 1)))
    ///     .collect();
    /// assert_eq!(
    ///     edges,
    ///     &[(0, 1), (1, 2), (2, 3), (3, 0),
    ///       (4, 5), (5, 6), (6, 7), (7, 4),
    ///       (0, 4), (1, 5), (2, 6), (3, 7),
    ///     ]
    /// );
    /// ```
    pub fn edge_node_id(&self, e: usize, i: usize) -> usize {
        match self {
            // Lin
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => 0,
            GeoKind::Lin4 => 0,
            GeoKind::Lin5 => 0,
            // Tri
            GeoKind::Tri3 => Tri3::EDGE_NODE_IDS[e][i],
            GeoKind::Tri6 => Tri6::EDGE_NODE_IDS[e][i],
            GeoKind::Tri10 => Tri10::EDGE_NODE_IDS[e][i],
            GeoKind::Tri15 => Tri15::EDGE_NODE_IDS[e][i],
            // Qua
            GeoKind::Qua4 => Qua4::EDGE_NODE_IDS[e][i],
            GeoKind::Qua8 => Qua8::EDGE_NODE_IDS[e][i],
            GeoKind::Qua9 => Qua9::EDGE_NODE_IDS[e][i],
            GeoKind::Qua12 => Qua12::EDGE_NODE_IDS[e][i],
            GeoKind::Qua16 => Qua16::EDGE_NODE_IDS[e][i],
            GeoKind::Qua17 => Qua17::EDGE_NODE_IDS[e][i],
            // Tet
            GeoKind::Tet4 => Tet4::EDGE_NODE_IDS[e][i],
            GeoKind::Tet10 => Tet10::EDGE_NODE_IDS[e][i],
            GeoKind::Tet20 => Tet20::EDGE_NODE_IDS[e][i],
            // Hex
            GeoKind::Hex8 => Hex8::EDGE_NODE_IDS[e][i],
            GeoKind::Hex20 => Hex20::EDGE_NODE_IDS[e][i],
            GeoKind::Hex32 => Hex32::EDGE_NODE_IDS[e][i],
        }
    }

    /// Returns the local id of node on edge (inward normals)
    ///
    /// # Input
    ///
    /// * `e` -- index of edge in [0, nedge-1]
    /// * `i` -- index of local node [0, edge_nnode-1]
    ///
    /// **Note:** For 2D GeoKinds, the sequence of edge nodes is such that the normals are **inward**.
    /// For 3D GeoKinds, the sequence just follows the associate Lin GeoKind.
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::GeoKind;
    ///
    /// //         2            OUTWARD NORMALS
    /// //   3-----6-----2          p0 p1 p2
    /// //   |           |      e:0 [1, 0, 4]
    /// //   |           |      e:1 [2, 1, 5]
    /// // 3 7           5 1    e:2 [3, 2, 6]
    /// //   |           |      e:3 [0, 3, 7]
    /// //   |           |
    /// //   0-----4-----1
    /// //         0
    /// let kind = GeoKind::Qua8;
    /// let edges: Vec<_> = (0..kind.nedge())
    ///     .map(|e| (kind.edge_node_id_inward(e, 0),
    ///               kind.edge_node_id_inward(e, 1),
    ///               kind.edge_node_id_inward(e, 2)))
    ///     .collect();
    /// assert_eq!(
    ///     edges,
    ///     &[(0, 1, 4), (1, 2, 5), (2, 3, 6), (3, 0, 7)]
    /// );
    /// ```
    ///
    /// ```
    /// use gemlab::shapes::GeoKind;
    ///
    /// //          .4-----[7]------7
    /// //        ,' |            ,'|
    /// //      [4] [8]         [6] |  [#] indicates local
    /// //    ,'     |        ,'    |      edge number "e"
    /// //  5'========[5]===6'     [11]
    /// //  |               |       |
    /// //  |        |      |       |
    /// // [9]      ,0-[3]- | - - - 3
    /// //  |     ,'       [10]  [2]
    /// //  |   [0]         |   ,'
    /// //  | ,'            | ,'
    /// //  1'-----[1]------2'
    /// let kind = GeoKind::Hex8;
    /// let edges: Vec<_> = (0..kind.nedge())
    ///     .map(|e| (kind.edge_node_id_inward(e, 0), kind.edge_node_id_inward(e, 1)))
    ///     .collect();
    /// assert_eq!(
    ///     edges,
    ///     &[(0, 1), (1, 2), (2, 3), (3, 0),
    ///       (4, 5), (5, 6), (6, 7), (7, 4),
    ///       (0, 4), (1, 5), (2, 6), (3, 7),
    ///     ]
    /// );
    /// ```
    pub fn edge_node_id_inward(&self, e: usize, i: usize) -> usize {
        match self {
            // Lin
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => 0,
            GeoKind::Lin4 => 0,
            GeoKind::Lin5 => 0,
            // Tri
            GeoKind::Tri3 => Tri3::EDGE_NODE_IDS_INWARD[e][i],
            GeoKind::Tri6 => Tri6::EDGE_NODE_IDS_INWARD[e][i],
            GeoKind::Tri10 => Tri10::EDGE_NODE_IDS_INWARD[e][i],
            GeoKind::Tri15 => Tri15::EDGE_NODE_IDS_INWARD[e][i],
            // Qua
            GeoKind::Qua4 => Qua4::EDGE_NODE_IDS_INWARD[e][i],
            GeoKind::Qua8 => Qua8::EDGE_NODE_IDS_INWARD[e][i],
            GeoKind::Qua9 => Qua9::EDGE_NODE_IDS_INWARD[e][i],
            GeoKind::Qua12 => Qua12::EDGE_NODE_IDS_INWARD[e][i],
            GeoKind::Qua16 => Qua16::EDGE_NODE_IDS_INWARD[e][i],
            GeoKind::Qua17 => Qua17::EDGE_NODE_IDS_INWARD[e][i],
            // Tet
            GeoKind::Tet4 => Tet4::EDGE_NODE_IDS[e][i],
            GeoKind::Tet10 => Tet10::EDGE_NODE_IDS[e][i],
            GeoKind::Tet20 => Tet20::EDGE_NODE_IDS[e][i],
            // Hex
            GeoKind::Hex8 => Hex8::EDGE_NODE_IDS[e][i],
            GeoKind::Hex20 => Hex20::EDGE_NODE_IDS[e][i],
            GeoKind::Hex32 => Hex32::EDGE_NODE_IDS[e][i],
        }
    }

    /// Returns the local id of node on face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `i` -- index of local node [0, face_nnode-1]
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::GeoKind;
    ///
    /// //           4----------------7
    /// //         ,'|              ,'|
    /// //       ,'  |  ___       ,'  |
    /// //     ,'    |,'5,'  [0],'    |
    /// //   ,'      |~~~     ,'      |
    /// // 5'===============6'  ,'|   |
    /// // |   ,'|   |      |   |3|   |
    /// // |   |2|   |      |   |,'   |
    /// // |   |,'   0- - - | +- - - -3
    /// // |       ,'       |       ,'
    /// // |     ,' [1]  ___|     ,'
    /// // |   ,'      ,'4,'|   ,'
    /// // | ,'        ~~~  | ,'
    /// // 1----------------2'
    /// let kind = GeoKind::Hex8;
    /// let faces: Vec<_> = (0..kind.nface())
    ///     .map(|f| {
    ///         (kind.face_node_id(f, 0), kind.face_node_id(f, 1),
    ///          kind.face_node_id(f, 2), kind.face_node_id(f, 3))
    ///     }).collect();
    /// assert_eq!(
    ///     faces,
    ///     &[
    ///         (0, 4, 7, 3), (1, 2, 6, 5),
    ///         (0, 1, 5, 4), (2, 3, 7, 6),
    ///         (0, 3, 2, 1), (4, 5, 6, 7),
    ///     ]
    /// );
    /// ```
    pub fn face_node_id(&self, f: usize, i: usize) -> usize {
        match self {
            // Lin
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => 0,
            GeoKind::Lin4 => 0,
            GeoKind::Lin5 => 0,
            // Tri
            GeoKind::Tri3 => 0,
            GeoKind::Tri6 => 0,
            GeoKind::Tri10 => 0,
            GeoKind::Tri15 => 0,
            // Qua
            GeoKind::Qua4 => 0,
            GeoKind::Qua8 => 0,
            GeoKind::Qua9 => 0,
            GeoKind::Qua12 => 0,
            GeoKind::Qua16 => 0,
            GeoKind::Qua17 => 0,
            // Tet
            GeoKind::Tet4 => Tet4::FACE_NODE_IDS[f][i],
            GeoKind::Tet10 => Tet10::FACE_NODE_IDS[f][i],
            GeoKind::Tet20 => Tet20::FACE_NODE_IDS[f][i],
            // Hex
            GeoKind::Hex8 => Hex8::FACE_NODE_IDS[f][i],
            GeoKind::Hex20 => Hex20::FACE_NODE_IDS[f][i],
            GeoKind::Hex32 => Hex32::FACE_NODE_IDS[f][i],
        }
    }

    /// Returns the local node id on an edge on the face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `k` -- index of face's edge (not the index of cell's edge) in [0, face_nedge-1]
    /// * `i` -- index of local node [0, edge_nnode-1]
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::GeoKind;
    ///
    /// //           4----------------7
    /// //         ,'|              ,'|
    /// //       ,'  |  ___       ,'  |
    /// //     ,'    |,'5,'  [0],'    |
    /// //   ,'      |~~~     ,'      |
    /// // 5'===============6'  ,'|   |
    /// // |   ,'|   |      |   |3|   |
    /// // |   |2|   |      |   |,'   |
    /// // |   |,'   0- - - | +- - - -3
    /// // |       ,'       |       ,'
    /// // |     ,' [1]  ___|     ,'
    /// // |   ,'      ,'4,'|   ,'
    /// // | ,'        ~~~  | ,'
    /// // 1----------------2'
    /// let kind = GeoKind::Hex8;
    /// let data: Vec<Vec<_>> = (0..kind.nface())
    ///     .map(|f| {
    ///         (0..kind.face_nedge())
    ///             .map(|k|
    ///                 (kind.face_edge_node_id(f, k, 0),
    ///                  kind.face_edge_node_id(f, k, 1))
    ///             ).collect()
    ///     }).collect();
    /// assert_eq!(
    ///     data,
    ///     &[
    ///         [(0, 4), (4, 7), (7, 3), (3, 0)], // face 0
    ///         [(1, 2), (2, 6), (6, 5), (5, 1)], // face 1
    ///         [(0, 1), (1, 5), (5, 4), (4, 0)], // face 2
    ///         [(2, 3), (3, 7), (7, 6), (6, 2)], // face 3
    ///         [(0, 3), (3, 2), (2, 1), (1, 0)], // face 4
    ///         [(4, 5), (5, 6), (6, 7), (7, 4)], // face 5
    ///     ]
    /// );
    /// ```
    pub fn face_edge_node_id(&self, f: usize, k: usize, i: usize) -> usize {
        match self {
            // Lin
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => 0,
            GeoKind::Lin4 => 0,
            GeoKind::Lin5 => 0,
            // Tri
            GeoKind::Tri3 => 0,
            GeoKind::Tri6 => 0,
            GeoKind::Tri10 => 0,
            GeoKind::Tri15 => 0,
            // Qua
            GeoKind::Qua4 => 0,
            GeoKind::Qua8 => 0,
            GeoKind::Qua9 => 0,
            GeoKind::Qua12 => 0,
            GeoKind::Qua16 => 0,
            GeoKind::Qua17 => 0,
            // Tet
            GeoKind::Tet4 => Tet4::FACE_EDGE_NODE_IDS[f][k][i],
            GeoKind::Tet10 => Tet10::FACE_EDGE_NODE_IDS[f][k][i],
            GeoKind::Tet20 => Tet20::FACE_EDGE_NODE_IDS[f][k][i],
            // Hex
            GeoKind::Hex8 => Hex8::FACE_EDGE_NODE_IDS[f][k][i],
            GeoKind::Hex20 => Hex20::FACE_EDGE_NODE_IDS[f][k][i],
            GeoKind::Hex32 => Hex32::FACE_EDGE_NODE_IDS[f][k][i],
        }
    }

    /// Returns the number of interior nodes (Lin{3,4,5}, Tri{10,15}, Qua{9,16,17})
    pub fn n_interior_nodes(&self) -> usize {
        match self {
            // Lin
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => Lin3::N_INTERIOR_NODE,
            GeoKind::Lin4 => Lin4::N_INTERIOR_NODE,
            GeoKind::Lin5 => Lin5::N_INTERIOR_NODE,
            // Tri
            GeoKind::Tri3 => 0,
            GeoKind::Tri6 => 0,
            GeoKind::Tri10 => Tri10::N_INTERIOR_NODE,
            GeoKind::Tri15 => Tri15::N_INTERIOR_NODE,
            // Qua
            GeoKind::Qua4 => 0,
            GeoKind::Qua8 => 0,
            GeoKind::Qua9 => Qua9::N_INTERIOR_NODE,
            GeoKind::Qua12 => 0,
            GeoKind::Qua16 => Qua16::N_INTERIOR_NODE,
            GeoKind::Qua17 => Qua17::N_INTERIOR_NODE,
            // Tet
            GeoKind::Tet4 => 0,
            GeoKind::Tet10 => 0,
            GeoKind::Tet20 => 0,
            // Hex
            GeoKind::Hex8 => 0,
            GeoKind::Hex20 => 0,
            GeoKind::Hex32 => 0,
        }
    }

    /// Returns the index of an interior node (Lin{3,4,5}, Tri{10,15}, Qua{9,16,17})
    ///
    /// **Note::** Make sure to call `n_interior_node` first to find out
    /// how many interior nodes a shape has; otherwise a panic may occur
    pub fn interior_node(&self, i: usize) -> usize {
        match self {
            // Lin
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => Lin3::INTERIOR_NODES[i],
            GeoKind::Lin4 => Lin4::INTERIOR_NODES[i],
            GeoKind::Lin5 => Lin5::INTERIOR_NODES[i],
            // Tri
            GeoKind::Tri3 => 0,
            GeoKind::Tri6 => 0,
            GeoKind::Tri10 => Tri10::INTERIOR_NODES[i],
            GeoKind::Tri15 => Tri15::INTERIOR_NODES[i],
            // Qua
            GeoKind::Qua4 => 0,
            GeoKind::Qua8 => 0,
            GeoKind::Qua9 => Qua9::INTERIOR_NODES[i],
            GeoKind::Qua12 => 0,
            GeoKind::Qua16 => Qua16::INTERIOR_NODES[i],
            GeoKind::Qua17 => Qua17::INTERIOR_NODES[i],
            // Tet
            GeoKind::Tet4 => 0,
            GeoKind::Tet10 => 0,
            GeoKind::Tet20 => 0,
            // Hex
            GeoKind::Hex8 => 0,
            GeoKind::Hex20 => 0,
            GeoKind::Hex32 => 0,
        }
    }

    /// Returns the reference coordinates at node m
    ///
    /// # Output
    ///
    /// * `ksi` -- (geo_ndim) reference coordinates `ξᵐ` at node m
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::GeoKind;
    ///
    /// //  3-------------2         ξ₀   ξ₁
    /// //  |      ξ₁     |  node    r    s
    /// //  |      |      |     0 -1.0 -1.0
    /// //  |      +--ξ₀  |     1  1.0 -1.0
    /// //  |             |     2  1.0  1.0
    /// //  |             |     3 -1.0  1.0
    /// //  0-------------1
    /// let kind = GeoKind::Qua4;
    /// let ref_coords: Vec<_> = (0..kind.nnode())
    ///     .map(|m| kind.reference_coords(m))
    ///     .collect();
    /// assert_eq!(ref_coords, &[
    ///     [-1.0, -1.0],
    ///     [ 1.0, -1.0],
    ///     [ 1.0,  1.0],
    ///     [-1.0,  1.0],
    /// ]);
    /// ```
    pub fn reference_coords(&self, m: usize) -> &'static [f64] {
        match self {
            // Lin
            GeoKind::Lin2 => &Lin2::NODE_REFERENCE_COORDS[m],
            GeoKind::Lin3 => &Lin3::NODE_REFERENCE_COORDS[m],
            GeoKind::Lin4 => &Lin4::NODE_REFERENCE_COORDS[m],
            GeoKind::Lin5 => &Lin5::NODE_REFERENCE_COORDS[m],
            // Tri
            GeoKind::Tri3 => &Tri3::NODE_REFERENCE_COORDS[m],
            GeoKind::Tri6 => &Tri6::NODE_REFERENCE_COORDS[m],
            GeoKind::Tri10 => &Tri10::NODE_REFERENCE_COORDS[m],
            GeoKind::Tri15 => &Tri15::NODE_REFERENCE_COORDS[m],
            // Qua
            GeoKind::Qua4 => &Qua4::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua8 => &Qua8::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua9 => &Qua9::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua12 => &Qua12::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua16 => &Qua16::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua17 => &Qua17::NODE_REFERENCE_COORDS[m],
            // Tet
            GeoKind::Tet4 => &Tet4::NODE_REFERENCE_COORDS[m],
            GeoKind::Tet10 => &Tet10::NODE_REFERENCE_COORDS[m],
            GeoKind::Tet20 => &Tet20::NODE_REFERENCE_COORDS[m],
            // Hex
            GeoKind::Hex8 => &Hex8::NODE_REFERENCE_COORDS[m],
            GeoKind::Hex20 => &Hex20::NODE_REFERENCE_COORDS[m],
            GeoKind::Hex32 => &Hex32::NODE_REFERENCE_COORDS[m],
        }
    }

    /// Returns the (min,delta) values on the reference domain
    ///
    /// # Tri and Tet
    ///
    /// * `ξ_min` = 0.0
    /// * `ξ_max` = 1.0
    /// * `Δξ` = 1.0
    ///
    /// # Lin, Qua, Hex
    ///
    /// * `ξ_min` = -1.0
    /// * `ξ_max` = +1.0
    /// * `Δξ` = 2.0
    pub fn ksi_min_ksi_del(&self) -> (f64, f64) {
        match self {
            // Lin
            GeoKind::Lin2 => (-1.0, 2.0),
            GeoKind::Lin3 => (-1.0, 2.0),
            GeoKind::Lin4 => (-1.0, 2.0),
            GeoKind::Lin5 => (-1.0, 2.0),
            // Tri
            GeoKind::Tri3 => (0.0, 1.0),
            GeoKind::Tri6 => (0.0, 1.0),
            GeoKind::Tri10 => (0.0, 1.0),
            GeoKind::Tri15 => (0.0, 1.0),
            // Qua
            GeoKind::Qua4 => (-1.0, 2.0),
            GeoKind::Qua8 => (-1.0, 2.0),
            GeoKind::Qua9 => (-1.0, 2.0),
            GeoKind::Qua12 => (-1.0, 2.0),
            GeoKind::Qua16 => (-1.0, 2.0),
            GeoKind::Qua17 => (-1.0, 2.0),
            // Tet
            GeoKind::Tet4 => (0.0, 1.0),
            GeoKind::Tet10 => (0.0, 1.0),
            GeoKind::Tet20 => (0.0, 1.0),
            // Hex
            GeoKind::Hex8 => (-1.0, 2.0),
            GeoKind::Hex20 => (-1.0, 2.0),
            GeoKind::Hex32 => (-1.0, 2.0),
        }
    }

    /// Returns the functions to evaluate the interpolation functions
    /// and the derivatives of the interpolation functions
    pub fn functions(&self) -> (FnInterp, FnDeriv) {
        match self {
            // Lin
            Self::Lin2 => (Lin2::calc_interp, Lin2::calc_deriv),
            Self::Lin3 => (Lin3::calc_interp, Lin3::calc_deriv),
            Self::Lin4 => (Lin4::calc_interp, Lin4::calc_deriv),
            Self::Lin5 => (Lin5::calc_interp, Lin5::calc_deriv),
            // Tri
            Self::Tri3 => (Tri3::calc_interp, Tri3::calc_deriv),
            Self::Tri6 => (Tri6::calc_interp, Tri6::calc_deriv),
            Self::Tri10 => (Tri10::calc_interp, Tri10::calc_deriv),
            Self::Tri15 => (Tri15::calc_interp, Tri15::calc_deriv),
            // Qua
            Self::Qua4 => (Qua4::calc_interp, Qua4::calc_deriv),
            Self::Qua8 => (Qua8::calc_interp, Qua8::calc_deriv),
            Self::Qua9 => (Qua9::calc_interp, Qua9::calc_deriv),
            Self::Qua12 => (Qua12::calc_interp, Qua12::calc_deriv),
            Self::Qua16 => (Qua16::calc_interp, Qua16::calc_deriv),
            Self::Qua17 => (Qua17::calc_interp, Qua17::calc_deriv),
            // Tet
            Self::Tet4 => (Tet4::calc_interp, Tet4::calc_deriv),
            Self::Tet10 => (Tet10::calc_interp, Tet10::calc_deriv),
            Self::Tet20 => (Tet20::calc_interp, Tet20::calc_deriv),
            // Hex
            Self::Hex8 => (Hex8::calc_interp, Hex8::calc_deriv),
            Self::Hex20 => (Hex20::calc_interp, Hex20::calc_deriv),
            Self::Hex32 => (Hex32::calc_interp, Hex32::calc_deriv),
        }
    }

    /// Returns true if Lin
    pub fn is_lin(&self) -> bool {
        match self {
            // Lin
            GeoKind::Lin2 => true,
            GeoKind::Lin3 => true,
            GeoKind::Lin4 => true,
            GeoKind::Lin5 => true,
            // Tri
            GeoKind::Tri3 => false,
            GeoKind::Tri6 => false,
            GeoKind::Tri10 => false,
            GeoKind::Tri15 => false,
            // Qua
            GeoKind::Qua4 => false,
            GeoKind::Qua8 => false,
            GeoKind::Qua9 => false,
            GeoKind::Qua12 => false,
            GeoKind::Qua16 => false,
            GeoKind::Qua17 => false,
            // Tet
            GeoKind::Tet4 => false,
            GeoKind::Tet10 => false,
            GeoKind::Tet20 => false,
            // Hex
            GeoKind::Hex8 => false,
            GeoKind::Hex20 => false,
            GeoKind::Hex32 => false,
        }
    }

    /// Returns true if Tri or Tet
    pub fn is_tri_or_tet(&self) -> bool {
        match self {
            // Lin
            GeoKind::Lin2 => false,
            GeoKind::Lin3 => false,
            GeoKind::Lin4 => false,
            GeoKind::Lin5 => false,
            // Tri
            GeoKind::Tri3 => true,
            GeoKind::Tri6 => true,
            GeoKind::Tri10 => true,
            GeoKind::Tri15 => true,
            // Qua
            GeoKind::Qua4 => false,
            GeoKind::Qua8 => false,
            GeoKind::Qua9 => false,
            GeoKind::Qua12 => false,
            GeoKind::Qua16 => false,
            GeoKind::Qua17 => false,
            // Tet
            GeoKind::Tet4 => true,
            GeoKind::Tet10 => true,
            GeoKind::Tet20 => true,
            // Hex
            GeoKind::Hex8 => false,
            GeoKind::Hex20 => false,
            GeoKind::Hex32 => false,
        }
    }

    /// Returns true if Qua or Hex
    pub fn is_qua_or_hex(&self) -> bool {
        match self {
            // Lin
            GeoKind::Lin2 => false,
            GeoKind::Lin3 => false,
            GeoKind::Lin4 => false,
            GeoKind::Lin5 => false,
            // Tri
            GeoKind::Tri3 => false,
            GeoKind::Tri6 => false,
            GeoKind::Tri10 => false,
            GeoKind::Tri15 => false,
            // Qua
            GeoKind::Qua4 => true,
            GeoKind::Qua8 => true,
            GeoKind::Qua9 => true,
            GeoKind::Qua12 => true,
            GeoKind::Qua16 => true,
            GeoKind::Qua17 => true,
            // Tet
            GeoKind::Tet4 => false,
            GeoKind::Tet10 => false,
            GeoKind::Tet20 => false,
            // Hex
            GeoKind::Hex8 => true,
            GeoKind::Hex20 => true,
            GeoKind::Hex32 => true,
        }
    }

    /// Returns the VTK cell type (if available in VTK)
    pub fn vtk_type(&self) -> Option<usize> {
        const VTK_LINE: usize = 3;
        const VTK_TRIANGLE: usize = 5;
        const VTK_QUAD: usize = 9;
        const VTK_TETRA: usize = 10;
        const VTK_HEXAHEDRON: usize = 12;
        const VTK_QUADRATIC_EDGE: usize = 21;
        const VTK_QUADRATIC_TRIANGLE: usize = 22;
        const VTK_QUADRATIC_QUAD: usize = 23;
        const VTK_QUADRATIC_TETRA: usize = 24;
        const VTK_QUADRATIC_HEXAHEDRON: usize = 25;
        const VTK_BIQUADRATIC_QUAD: usize = 28;
        const VTK_CUBIC_LINE: usize = 35;
        match self {
            // Lin
            GeoKind::Lin2 => Some(VTK_LINE),
            GeoKind::Lin3 => Some(VTK_QUADRATIC_EDGE),
            GeoKind::Lin4 => Some(VTK_CUBIC_LINE),
            GeoKind::Lin5 => None,
            // Tri
            GeoKind::Tri3 => Some(VTK_TRIANGLE),
            GeoKind::Tri6 => Some(VTK_QUADRATIC_TRIANGLE),
            GeoKind::Tri10 => None,
            GeoKind::Tri15 => None,
            // Qua
            GeoKind::Qua4 => Some(VTK_QUAD),
            GeoKind::Qua8 => Some(VTK_QUADRATIC_QUAD),
            GeoKind::Qua9 => Some(VTK_BIQUADRATIC_QUAD),
            GeoKind::Qua12 => None,
            GeoKind::Qua16 => None,
            GeoKind::Qua17 => None,
            // Tet
            GeoKind::Tet4 => Some(VTK_TETRA),
            GeoKind::Tet10 => Some(VTK_QUADRATIC_TETRA),
            GeoKind::Tet20 => None,
            // Hex
            GeoKind::Hex8 => Some(VTK_HEXAHEDRON),
            GeoKind::Hex20 => Some(VTK_QUADRATIC_HEXAHEDRON),
            GeoKind::Hex32 => None,
        }
    }

    /// Holds all enum values
    pub const VALUES: [Self; 20] = [
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
        Self::Tet20,
        // Hex
        Self::Hex8,
        Self::Hex20,
        Self::Hex32,
    ];
}

/// Converts edge-local-index to face-local-indices of a Hexahedron
pub const HEX_EDGE_TO_FACE: [[usize; 2]; 12] = [
    [2, 4], //  0
    [1, 4], //  1
    [3, 4], //  2
    [0, 4], //  3
    [2, 5], //  4
    [1, 5], //  5
    [3, 5], //  6
    [0, 5], //  7
    [0, 2], //  8
    [1, 2], //  9
    [1, 3], // 10
    [0, 3], // 11
];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{GeoCase, GeoClass, GeoKind};
    use std::collections::HashSet;

    const IGNORED: usize = 0;
    const SWAP_EDGE_NODE_LIN2: [usize; 5] = [1, 0, IGNORED, IGNORED, IGNORED];
    const SWAP_EDGE_NODE_LIN3: [usize; 5] = [1, 0, 2, IGNORED, IGNORED];
    const SWAP_EDGE_NODE_LIN4: [usize; 5] = [1, 0, 3, 2, IGNORED];
    const SWAP_EDGE_NODE_LIN5: [usize; 5] = [1, 0, 2, 4, 3];

    fn swap_edge_node_id(tri_or_qua_kind: GeoKind, e: usize, i: usize) -> usize {
        let edge_kind = tri_or_qua_kind.edge_kind().unwrap();
        let swap = match edge_kind {
            GeoKind::Lin2 => SWAP_EDGE_NODE_LIN2,
            GeoKind::Lin3 => SWAP_EDGE_NODE_LIN3,
            GeoKind::Lin4 => SWAP_EDGE_NODE_LIN4,
            GeoKind::Lin5 => SWAP_EDGE_NODE_LIN5,
            _ => panic!("INTERNAL ERROR: edge kind is not possible for Tri or Qua"),
        };
        tri_or_qua_kind.edge_node_id(e, swap[i])
    }

    #[test]
    fn derive_works() {
        let case = GeoCase::Cable.clone();
        let class = GeoClass::Tri.clone();
        let kind = GeoKind::Tri6.clone();
        let classes = HashSet::from([GeoClass::Tri, GeoClass::Qua]);
        let kinds = HashSet::from([GeoKind::Tri3, GeoKind::Qua4]);
        assert_eq!(case, GeoCase::Cable);
        assert_eq!(class, GeoClass::Tri);
        assert_eq!(kind, GeoKind::Tri6);
        assert_eq!(format!("{:?}", case), "Cable");
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
                    assert_eq!(GeoKind::from("lin2").unwrap(), kind);
                    assert_eq!(kind.to_string(), "lin2");
                    assert_eq!(kind.class(), GeoClass::Lin);
                    assert_eq!(kind.ndim(), 1);
                    assert_eq!(kind.nnode(), 2);
                    assert_eq!(kind.nedge(), 0);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 0);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), None);
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), Some(3));
                }
                GeoKind::Lin3 => {
                    assert_eq!(GeoKind::from("lin3").unwrap(), kind);
                    assert_eq!(kind.to_string(), "lin3");
                    assert_eq!(kind.class(), GeoClass::Lin);
                    assert_eq!(kind.ndim(), 1);
                    assert_eq!(kind.nnode(), 3);
                    assert_eq!(kind.nedge(), 0);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 0);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), None);
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), Some(GeoKind::Lin2));
                    assert_eq!(kind.n_interior_nodes(), 1);
                    assert_eq!(kind.interior_node(0), 2);
                    assert_eq!(kind.vtk_type(), Some(21));
                }
                GeoKind::Lin4 => {
                    assert_eq!(GeoKind::from("lin4").unwrap(), kind);
                    assert_eq!(kind.to_string(), "lin4");
                    assert_eq!(kind.class(), GeoClass::Lin);
                    assert_eq!(kind.ndim(), 1);
                    assert_eq!(kind.nnode(), 4);
                    assert_eq!(kind.nedge(), 0);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 0);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), None);
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 2);
                    assert_eq!(kind.interior_node(0), 2);
                    assert_eq!(kind.interior_node(1), 3);
                    assert_eq!(kind.vtk_type(), Some(35));
                }
                GeoKind::Lin5 => {
                    assert_eq!(GeoKind::from("lin5").unwrap(), kind);
                    assert_eq!(kind.to_string(), "lin5");
                    assert_eq!(kind.class(), GeoClass::Lin);
                    assert_eq!(kind.ndim(), 1);
                    assert_eq!(kind.nnode(), 5);
                    assert_eq!(kind.nedge(), 0);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 0);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), None);
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), Some(GeoKind::Lin3));
                    assert_eq!(kind.n_interior_nodes(), 3);
                    assert_eq!(kind.interior_node(0), 2);
                    assert_eq!(kind.interior_node(1), 3);
                    assert_eq!(kind.interior_node(2), 4);
                    assert_eq!(kind.vtk_type(), None);
                }

                // Tri
                GeoKind::Tri3 => {
                    assert_eq!(GeoKind::from("tri3").unwrap(), kind);
                    assert_eq!(kind.to_string(), "tri3");
                    assert_eq!(kind.class(), GeoClass::Tri);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 3);
                    assert_eq!(kind.nedge(), 3);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 2);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin2));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), Some(5));
                }
                GeoKind::Tri6 => {
                    assert_eq!(GeoKind::from("tri6").unwrap(), kind);
                    assert_eq!(kind.to_string(), "tri6");
                    assert_eq!(kind.class(), GeoClass::Tri);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 6);
                    assert_eq!(kind.nedge(), 3);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 3);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin3));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), Some(GeoKind::Tri3));
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), Some(22));
                }
                GeoKind::Tri10 => {
                    assert_eq!(GeoKind::from("tri10").unwrap(), kind);
                    assert_eq!(kind.to_string(), "tri10");
                    assert_eq!(kind.class(), GeoClass::Tri);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 10);
                    assert_eq!(kind.nedge(), 3);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 4);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin4));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 1);
                    assert_eq!(kind.interior_node(0), 9);
                    assert_eq!(kind.vtk_type(), None);
                }
                GeoKind::Tri15 => {
                    assert_eq!(GeoKind::from("tri15").unwrap(), kind);
                    assert_eq!(kind.to_string(), "tri15");
                    assert_eq!(kind.class(), GeoClass::Tri);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 15);
                    assert_eq!(kind.nedge(), 3);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 5);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin5));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), Some(GeoKind::Tri6));
                    assert_eq!(kind.n_interior_nodes(), 3);
                    assert_eq!(kind.interior_node(0), 12);
                    assert_eq!(kind.interior_node(1), 13);
                    assert_eq!(kind.interior_node(2), 14);
                    assert_eq!(kind.vtk_type(), None);
                }

                // Qua
                GeoKind::Qua4 => {
                    assert_eq!(GeoKind::from("qua4").unwrap(), kind);
                    assert_eq!(kind.to_string(), "qua4");
                    assert_eq!(kind.class(), GeoClass::Qua);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 4);
                    assert_eq!(kind.nedge(), 4);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 2);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin2));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), Some(9));
                }
                GeoKind::Qua8 => {
                    assert_eq!(GeoKind::from("qua8").unwrap(), kind);
                    assert_eq!(kind.to_string(), "qua8");
                    assert_eq!(kind.class(), GeoClass::Qua);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 8);
                    assert_eq!(kind.nedge(), 4);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 3);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin3));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), Some(GeoKind::Qua4));
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), Some(23));
                }
                GeoKind::Qua9 => {
                    assert_eq!(GeoKind::from("qua9").unwrap(), kind);
                    assert_eq!(kind.to_string(), "qua9");
                    assert_eq!(kind.class(), GeoClass::Qua);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 9);
                    assert_eq!(kind.nedge(), 4);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 3);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin3));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), Some(GeoKind::Qua4));
                    assert_eq!(kind.n_interior_nodes(), 1);
                    assert_eq!(kind.interior_node(0), 8);
                    assert_eq!(kind.vtk_type(), Some(28));
                }
                GeoKind::Qua12 => {
                    assert_eq!(GeoKind::from("qua12").unwrap(), kind);
                    assert_eq!(kind.to_string(), "qua12");
                    assert_eq!(kind.class(), GeoClass::Qua);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 12);
                    assert_eq!(kind.nedge(), 4);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 4);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin4));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), None);
                }
                GeoKind::Qua16 => {
                    assert_eq!(GeoKind::from("qua16").unwrap(), kind);
                    assert_eq!(kind.to_string(), "qua16");
                    assert_eq!(kind.class(), GeoClass::Qua);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 16);
                    assert_eq!(kind.nedge(), 4);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 4);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin4));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 4);
                    assert_eq!(kind.interior_node(0), 12);
                    assert_eq!(kind.interior_node(1), 13);
                    assert_eq!(kind.interior_node(2), 14);
                    assert_eq!(kind.interior_node(3), 15);
                    assert_eq!(kind.vtk_type(), None);
                }
                GeoKind::Qua17 => {
                    assert_eq!(GeoKind::from("qua17").unwrap(), kind);
                    assert_eq!(kind.to_string(), "qua17");
                    assert_eq!(kind.class(), GeoClass::Qua);
                    assert_eq!(kind.ndim(), 2);
                    assert_eq!(kind.nnode(), 17);
                    assert_eq!(kind.nedge(), 4);
                    assert_eq!(kind.nface(), 0);
                    assert_eq!(kind.edge_nnode(), 5);
                    assert_eq!(kind.face_nnode(), 0);
                    assert_eq!(kind.face_nedge(), 0);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin5));
                    assert_eq!(kind.face_kind(), None);
                    assert_eq!(kind.lower_order(), Some(GeoKind::Qua9));
                    assert_eq!(kind.n_interior_nodes(), 1);
                    assert_eq!(kind.interior_node(0), 8);
                    assert_eq!(kind.vtk_type(), None);
                }

                // Tet
                GeoKind::Tet4 => {
                    assert_eq!(GeoKind::from("tet4").unwrap(), kind);
                    assert_eq!(kind.to_string(), "tet4");
                    assert_eq!(kind.class(), GeoClass::Tet);
                    assert_eq!(kind.ndim(), 3);
                    assert_eq!(kind.nnode(), 4);
                    assert_eq!(kind.nedge(), 6);
                    assert_eq!(kind.nface(), 4);
                    assert_eq!(kind.edge_nnode(), 2);
                    assert_eq!(kind.face_nnode(), 3);
                    assert_eq!(kind.face_nedge(), 3);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin2));
                    assert_eq!(kind.face_kind(), Some(GeoKind::Tri3));
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), Some(10));
                }
                GeoKind::Tet10 => {
                    assert_eq!(GeoKind::from("tet10").unwrap(), kind);
                    assert_eq!(kind.to_string(), "tet10");
                    assert_eq!(kind.class(), GeoClass::Tet);
                    assert_eq!(kind.ndim(), 3);
                    assert_eq!(kind.nnode(), 10);
                    assert_eq!(kind.nedge(), 6);
                    assert_eq!(kind.nface(), 4);
                    assert_eq!(kind.edge_nnode(), 3);
                    assert_eq!(kind.face_nnode(), 6);
                    assert_eq!(kind.face_nedge(), 3);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin3));
                    assert_eq!(kind.face_kind(), Some(GeoKind::Tri6));
                    assert_eq!(kind.lower_order(), Some(GeoKind::Tet4));
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), Some(24));
                }
                GeoKind::Tet20 => {
                    assert_eq!(GeoKind::from("tet20").unwrap(), kind);
                    assert_eq!(kind.to_string(), "tet20");
                    assert_eq!(kind.class(), GeoClass::Tet);
                    assert_eq!(kind.ndim(), 3);
                    assert_eq!(kind.nnode(), 20);
                    assert_eq!(kind.nedge(), 6);
                    assert_eq!(kind.nface(), 4);
                    assert_eq!(kind.edge_nnode(), 4);
                    assert_eq!(kind.face_nnode(), 10);
                    assert_eq!(kind.face_nedge(), 3);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin4));
                    assert_eq!(kind.face_kind(), Some(GeoKind::Tri10));
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), None);
                }

                // Hex
                GeoKind::Hex8 => {
                    assert_eq!(GeoKind::from("hex8").unwrap(), kind);
                    assert_eq!(kind.to_string(), "hex8");
                    assert_eq!(kind.class(), GeoClass::Hex);
                    assert_eq!(kind.ndim(), 3);
                    assert_eq!(kind.nnode(), 8);
                    assert_eq!(kind.nedge(), 12);
                    assert_eq!(kind.nface(), 6);
                    assert_eq!(kind.edge_nnode(), 2);
                    assert_eq!(kind.face_nnode(), 4);
                    assert_eq!(kind.face_nedge(), 4);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin2));
                    assert_eq!(kind.face_kind(), Some(GeoKind::Qua4));
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), Some(12));
                }
                GeoKind::Hex20 => {
                    assert_eq!(GeoKind::from("hex20").unwrap(), kind);
                    assert_eq!(kind.to_string(), "hex20");
                    assert_eq!(kind.class(), GeoClass::Hex);
                    assert_eq!(kind.ndim(), 3);
                    assert_eq!(kind.nnode(), 20);
                    assert_eq!(kind.nedge(), 12);
                    assert_eq!(kind.nface(), 6);
                    assert_eq!(kind.edge_nnode(), 3);
                    assert_eq!(kind.face_nnode(), 8);
                    assert_eq!(kind.face_nedge(), 4);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin3));
                    assert_eq!(kind.face_kind(), Some(GeoKind::Qua8));
                    assert_eq!(kind.lower_order(), Some(GeoKind::Hex8));
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), Some(25));
                }
                GeoKind::Hex32 => {
                    assert_eq!(GeoKind::from("hex32").unwrap(), kind);
                    assert_eq!(kind.to_string(), "hex32");
                    assert_eq!(kind.class(), GeoClass::Hex);
                    assert_eq!(kind.ndim(), 3);
                    assert_eq!(kind.nnode(), 32);
                    assert_eq!(kind.nedge(), 12);
                    assert_eq!(kind.nface(), 6);
                    assert_eq!(kind.edge_nnode(), 4);
                    assert_eq!(kind.face_nnode(), 12);
                    assert_eq!(kind.face_nedge(), 4);
                    assert_eq!(kind.edge_kind(), Some(GeoKind::Lin4));
                    assert_eq!(kind.face_kind(), Some(GeoKind::Qua12));
                    assert_eq!(kind.lower_order(), None);
                    assert_eq!(kind.n_interior_nodes(), 0);
                    assert_eq!(kind.interior_node(0), 0);
                    assert_eq!(kind.vtk_type(), None);
                }
            }
            match kind.class() {
                GeoClass::Lin => {
                    assert_eq!(kind.edge_node_id(0, 0), 0);
                    assert_eq!(kind.edge_node_id_inward(0, 0), 0);
                    assert_eq!(kind.is_lin(), true);
                    assert_eq!(kind.is_tri_or_tet(), false);
                    assert_eq!(kind.is_qua_or_hex(), false);
                }
                GeoClass::Tri => {
                    assert_eq!(kind.edge_node_id(0, 0), 1); // it always going to be the next node, after 0
                    assert_eq!(kind.edge_node_id_inward(0, 0), 0);
                    assert_eq!(kind.is_lin(), false);
                    assert_eq!(kind.is_tri_or_tet(), true);
                    assert_eq!(kind.is_qua_or_hex(), false);
                    for e in 0..3 {
                        let inward: Vec<_> = (0..kind.edge_nnode()).map(|i| kind.edge_node_id_inward(e, i)).collect();
                        let correct: Vec<_> = (0..kind.edge_nnode()).map(|i| swap_edge_node_id(kind, e, i)).collect();
                        assert_eq!(inward, correct);
                    }
                }
                GeoClass::Qua => {
                    assert_eq!(kind.edge_node_id(0, 0), 1); // it always going to be the next node, after 0
                    assert_eq!(kind.edge_node_id_inward(0, 0), 0);
                    assert_eq!(kind.is_lin(), false);
                    assert_eq!(kind.is_tri_or_tet(), false);
                    assert_eq!(kind.is_qua_or_hex(), true);
                    for e in 0..4 {
                        let inward: Vec<_> = (0..kind.edge_nnode()).map(|i| kind.edge_node_id_inward(e, i)).collect();
                        let correct: Vec<_> = (0..kind.edge_nnode()).map(|i| swap_edge_node_id(kind, e, i)).collect();
                        assert_eq!(inward, correct);
                    }
                }
                GeoClass::Tet => {
                    assert_eq!(kind.edge_node_id(0, 1), 1); // the next node is equal to 1
                    assert_eq!(kind.edge_node_id_inward(0, 0), 0);
                    assert_eq!(kind.is_lin(), false);
                    assert_eq!(kind.is_tri_or_tet(), true);
                    assert_eq!(kind.is_qua_or_hex(), false);
                }
                GeoClass::Hex => {
                    assert_eq!(kind.edge_node_id(0, 1), 1); // the next node is equal to 1
                    assert_eq!(kind.edge_node_id_inward(0, 0), 0);
                    assert_eq!(kind.is_lin(), false);
                    assert_eq!(kind.is_tri_or_tet(), false);
                    assert_eq!(kind.is_qua_or_hex(), true);
                }
            }
            assert_eq!(kind.face_node_id(0, 0), 0);
            assert_eq!(kind.face_edge_node_id(0, 0, 0), 0);
            if let Some(lo) = kind.lower_order() {
                for m in 0..lo.nnode() {
                    assert_eq!(lo.reference_coords(m), kind.reference_coords(m));
                }
            }
        }
    }
}

use super::{AsCell, CellId, PointId};
use crate::shapes::GeoKind;
use russell_lab::sort4;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Defines a unique key for faces by using a quadruple of sorted point indices
///
/// **Note:** For 3-node faces, the fourth entry in the key will be set to [usize::MAX].
/// In this way, we can compare 4-node (or more nodes) faces. Since the local numbering
/// scheme runs over the "corners" first, the middle points don't matter.
pub type FaceKey = (usize, usize, usize, usize);

/// Holds the essential information to reconstruct an face
///
/// * A face is an entity belonging to a solid cell in 3D
#[derive(Clone, Debug)]
pub struct Face {
    /// Geometry kind
    pub kind: GeoKind,

    /// List of points defining this edge or face; in the right (FEM) order (i.e., unsorted)
    pub points: Vec<PointId>,

    /// Marker
    pub marker: i32,
}

impl Face {
    /// Returns the sorted list of key points
    pub fn key(&self) -> FaceKey {
        let mut key = if self.points.len() > 3 {
            (self.points[0], self.points[1], self.points[2], self.points[3])
        } else {
            (self.points[0], self.points[1], self.points[2], usize::MAX)
        };
        sort4(&mut key);
        key
    }
}

/// Defines an array of faces
#[derive(Clone, Debug)]
pub struct Faces<'a> {
    /// Holds a set of faces
    pub all: Vec<&'a Face>,
}

/// Maps faces to cells sharing the face (3D only)
///
/// Relates face keys to `Vec<(cell_id, f)>` where:
///
/// * `cell_id` -- is the id of the cell sharing the face
/// * `f` -- is the cell's local face index
pub type MapFaceToCells = HashMap<FaceKey, Vec<(CellId, usize)>>;

/// Maps a point id to faces sharing the point
///
/// Relates a point id to a unique set of FaceKey
pub type MapPointToFaces = HashMap<PointId, HashSet<FaceKey>>;

impl AsCell for Face {
    fn kind(&self) -> GeoKind {
        self.kind
    }

    fn marker(&self) -> i32 {
        self.marker
    }

    fn points(&self) -> &[PointId] {
        &self.points
    }
}

impl fmt::Display for Face {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.points.len() == 3 {
            let (a, b, c, _) = self.key();
            write!(f, "({}, {}, {}, MAX)", a, b, c).unwrap();
        } else {
            write!(f, "{:?}", self.key()).unwrap();
        }
        Ok(())
    }
}

impl<'a> fmt::Display for Faces<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.all.len() {
            if i > 0 {
                write!(f, ", ").unwrap();
            }
            write!(f, "{}", self.all[i]).unwrap();
        }
        Ok(())
    }
}

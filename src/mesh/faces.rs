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

impl<'a> Faces<'a> {
    /// Returns a sorted list of all unique points in the collection of faces
    ///
    /// # Returns
    ///
    /// A vector containing all unique point IDs from all faces in `self.all`, sorted in ascending order.
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::{Face, Faces};
    /// use gemlab::shapes::GeoKind;
    ///
    /// let f1 = Face { kind: GeoKind::Tri3, points: vec![1, 2, 3], marker: 0 };
    /// let f2 = Face { kind: GeoKind::Tri3, points: vec![2, 3, 4], marker: 0 };
    /// let faces = Faces { all: vec![&f1, &f2] };
    ///
    /// let points = faces.all_points();
    /// assert_eq!(points, vec![1, 2, 3, 4]);
    /// ```
    pub fn all_points(&self) -> Vec<PointId> {
        let mut points_set = HashSet::new();
        for face in &self.all {
            for &point in &face.points {
                points_set.insert(point);
            }
        }
        let mut points: Vec<_> = points_set.into_iter().collect();
        points.sort();
        points
    }
}

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Face, Faces};
    use crate::shapes::GeoKind;

    #[rustfmt::skip]
    fn generate_sample_1() -> Vec<Face> {
        vec![
            Face { kind: GeoKind::Tri3, points: vec![1, 2, 3], marker: 0 },
            Face { kind: GeoKind::Tri3, points: vec![4, 5, 6], marker: 0 },
        ]
    }

    #[rustfmt::skip]
    fn generate_sample_2() -> Vec<Face> {
        vec![
            Face { kind: GeoKind::Qua4, points: vec![1, 2, 3, 4], marker: 0 },
            Face { kind: GeoKind::Qua4, points: vec![2, 5, 6, 3], marker: 0 },
            Face { kind: GeoKind::Tri6, points: vec![7, 8, 9, 10, 11, 12], marker: 0 },
        ]
    }

    #[test]
    fn all_points_works() {
        // Empty list of faces
        let empty = Faces { all: vec![] };
        assert!(empty.all_points().is_empty());

        // Single face (Tri3)
        let all = generate_sample_1();
        let single = Faces { all: vec![&all[0]] };
        assert_eq!(single.all_points(), vec![1, 2, 3]);

        // Two faces with no shared points
        let two_faces = Faces {
            all: vec![&all[0], &all[1]],
        };
        assert_eq!(two_faces.all_points(), vec![1, 2, 3, 4, 5, 6]);

        // Mixed face types with shared points
        let all2 = generate_sample_2();
        let mixed = Faces {
            all: vec![&all2[0], &all2[1]],
        };
        assert_eq!(mixed.all_points(), vec![1, 2, 3, 4, 5, 6]);

        // Higher-order face (Tri6 with middle nodes)
        let higher_order = Faces { all: vec![&all2[2]] };
        assert_eq!(higher_order.all_points(), vec![7, 8, 9, 10, 11, 12]);

        // All faces together
        let all_faces = Faces {
            all: vec![&all2[0], &all2[1], &all2[2]],
        };
        assert_eq!(all_faces.all_points(), vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
    }
}

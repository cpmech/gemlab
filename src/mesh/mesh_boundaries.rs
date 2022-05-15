use super::{alloc_cell_shapes, Cell, CellId, Mesh, PointId};
use crate::{shapes::Shape, StrError};
use russell_lab::{sort2, sort4};
use std::collections::{HashMap, HashSet};

/// Aliases (usize,usize) as the key of Edge
///
/// # Note
///
/// Since the local numbering scheme runs over "corners" first, we can compare
/// edges using only two points; i.e., the middle points don't matter.
pub type EdgeKey = (usize, usize);

/// Aliases (usize,usize,usize,usize) as the key of Face
///
/// # Note
///
/// If all faces have at most 3 points, the fourth entry in the key will be equal to the total number of points.
/// In this way, we can compare 4-node (or more nodes) faces with each other, since that the local numbering
/// scheme runs over the "corners" first; i.e., the middle points don't matter.
pub type FaceKey = (usize, usize, usize, usize);

/// Holds edge data (derived data structure)
#[derive(Clone, Debug)]
pub struct Edge {
    /// List of points defining this edge; in the right order (unsorted)
    pub points: Vec<PointId>,
}

/// Holds face data (derived data structure)
#[derive(Clone, Debug)]
pub struct Face {
    /// List of points defining this face; in the right order (unsorted)
    pub points: Vec<PointId>,
}

/// Holds points, edges and faces on the boundaries of a mesh
pub struct MeshBoundaries {
    /// Set of points on the boundaries
    ///
    /// Note: a boundary point belongs to a boundary edge or a boundary face
    pub points: HashSet<PointId>,

    /// Set of edges on the boundaries
    ///
    /// Note:
    ///
    /// * In 2D, a boundary edge is such that it is shared by one 2D cell only (1D cells are ignored)
    /// * In 3D, a boundary edge belongs to a boundary face
    pub edges: HashMap<EdgeKey, Edge>,

    /// Set of faces on the boundaries
    ///
    /// Note: A boundary face is such that it is shared by one 3D cell only
    pub faces: HashMap<FaceKey, Face>,
}

impl MeshBoundaries {
    /// Finds the boundary entities and allocates a new instance
    pub fn new(mesh: &Mesh) -> Result<Self, StrError> {
        // maps all face keys to (cell_id, f) where f is the cell's local face index
        // let mut all_faces: HashMap<FaceKey, Vec<(CellId, usize)>> = HashMap::new(); // (face_key) => [(cell_id,f)]

        // let shapes = alloc_cell_shapes(mesh)?;
        // let mut edges: HashMap<EdgeKey, Edge> = HashMap::new();
        // let mut faces: HashMap<FaceKey, Face> = HashMap::new();

        // match mesh.space_ndim {
        //     2 => find_boundaries_2d(mesh, &shapes)?,
        //     _ => panic!("space_ndim must be 2 or 3"),
        // };

        Ok(MeshBoundaries {
            points: HashSet::new(),
            edges: HashMap::new(),
            faces: HashMap::new(),
        })
    }
}

/// Maps all edges (internal and boundary) in 2D
///
/// # Output
///
/// Returns a map relating edge keys to `Vec<(cell_id, e)>` where:
///
/// * `cell_id` -- the id of the cell sharing the edge
/// * `e` -- is the cell's local edge index
///
/// # Panics
///
/// The `mesh.space_ndim` must be equal to 2, otherwise a panic occurs.
pub fn map_edges_2d(mesh: &Mesh, shapes: &Vec<Shape>) -> HashMap<EdgeKey, Vec<(CellId, usize)>> {
    assert_eq!(mesh.space_ndim, 2);
    let mut edges: HashMap<EdgeKey, Vec<(CellId, usize)>> = HashMap::new();
    mesh.cells.iter().zip(shapes).for_each(|(cell, shape)| {
        if shape.geo_ndim == 2 {
            for e in 0..shape.nedge {
                let mut edge_key: EdgeKey = (
                    cell.points[shape.edge_node_id(e, 0)],
                    cell.points[shape.edge_node_id(e, 1)],
                );
                sort2(&mut edge_key);
                let data = edges.entry(edge_key).or_insert(Vec::new());
                data.push((cell.id, e));
            }
        }
    });
    edges
}

/// Maps all faces (internal and boundary)
///
/// Returns a map relating face keys to `Vec<(cell_id, f)>` where:
///
/// * `cell_id` -- the id of the cell sharing the face
/// * `f` -- is the cell's local face index
///
/// # Panics
///
/// The `mesh.space_ndim` must be equal to 3, otherwise a panic occurs.
pub fn map_faces(mesh: &Mesh, shapes: &Vec<Shape>) -> HashMap<FaceKey, Vec<(CellId, usize)>> {
    assert_eq!(mesh.space_ndim, 3);
    let mut faces: HashMap<FaceKey, Vec<(CellId, usize)>> = HashMap::new();
    mesh.cells.iter().zip(shapes).for_each(|(cell, shape)| {
        if shape.geo_ndim == 3 {
            for f in 0..shape.nface {
                let mut face_key: FaceKey = if shape.face_nnode > 3 {
                    (
                        cell.points[shape.face_node_id(f, 0)],
                        cell.points[shape.face_node_id(f, 1)],
                        cell.points[shape.face_node_id(f, 2)],
                        cell.points[shape.face_node_id(f, 3)],
                    )
                } else {
                    (
                        cell.points[shape.face_node_id(f, 0)],
                        cell.points[shape.face_node_id(f, 1)],
                        cell.points[shape.face_node_id(f, 2)],
                        mesh.points.len(),
                    )
                };
                sort4(&mut face_key);
                let data = faces.entry(face_key).or_insert(Vec::new());
                data.push((cell.id, f));
            }
        }
    });
    faces
}

/*
fn find_boundaries_2d(
    mesh: &Mesh,
    shapes: &Vec<Shape>,
) -> Result<(HashSet<PointId>, HashMap<EdgeKey, Edge>), StrError> {
    // holds the ids of points on boundary or on shapes with geo_ndim == 1
    let mut points: HashSet<PointId> = HashSet::new();

    // find all edges
    mesh.cells
        .iter()
        .zip(shapes)
        .for_each(|(cell, shape)| match shape.geo_ndim {
            1 => {
                for m in 0..shape.nnode {
                    points.insert(cell.points[m]);
                }
            }
            2 => {
                for e in 0..shape.nedge {
                    let mut edge_key: EdgeKey = (
                        cell.points[shape.edge_node_id(e, 0)],
                        cell.points[shape.edge_node_id(e, 1)],
                    );
                    sort2(&mut edge_key);
                    match all_edges.get_mut(&edge_key) {
                        Some((_, _, internal)) => {
                            *internal = true;
                        }
                        None => {
                            all_edges.insert(edge_key, (cell.id, e, false));
                        }
                    }
                }
            }
            _ => panic!("shape.geo_ndim is invalid in 2D"),
        });

    // allocate boundary edges
    let edges: HashMap<_, _> = all_edges
        .iter()
        .filter_map(|(edge_key, (cell_id, e, internal))| match internal {
            true => None,
            false => Some((
                *edge_key,
                Edge {
                    points: (0..shapes[*cell_id].edge_nnode)
                        .into_iter()
                        .map(|i| shapes[*cell_id].edge_node_id(*e, i))
                        .collect(),
                },
            )),
        })
        .collect();

    Ok((points, edges))
}
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{map_edges_2d, map_faces};
    use crate::mesh::{alloc_cell_shapes, Samples};
    use crate::StrError;

    #[test]
    fn map_edges_2d_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let shapes = alloc_cell_shapes(&mesh)?;
        let edges = map_edges_2d(&mesh, &shapes);
        let mut keys: Vec<_> = edges.keys().collect();
        keys.sort();
        assert_eq!(keys, [&(0, 1), &(0, 3), &(1, 2), &(1, 4), &(2, 3), &(2, 5), &(4, 5)]);
        assert_eq!(edges.get(&(0, 1)).unwrap(), &[(0, 0)]);
        assert_eq!(edges.get(&(0, 3)).unwrap(), &[(0, 3)]);
        assert_eq!(edges.get(&(1, 2)).unwrap(), &[(0, 1), (1, 3)]);
        assert_eq!(edges.get(&(1, 4)).unwrap(), &[(1, 0)]);
        assert_eq!(edges.get(&(2, 3)).unwrap(), &[(0, 2)]);
        assert_eq!(edges.get(&(2, 5)).unwrap(), &[(1, 2)]);
        assert_eq!(edges.get(&(4, 5)).unwrap(), &[(1, 1)]);

        //           4---------3
        //           |         |
        //           |   [1]   |
        //           |         |
        //  0--------1---------2
        let mesh = Samples::mixed_shapes_2d();
        let shapes = alloc_cell_shapes(&mesh)?;
        let edges = map_edges_2d(&mesh, &shapes);
        let mut keys: Vec<_> = edges.keys().collect();
        keys.sort();
        assert_eq!(keys, [&(1, 2), &(1, 4), &(2, 3), &(3, 4)]);
        assert_eq!(edges.get(&(1, 2)).unwrap(), &[(1, 0)]);
        assert_eq!(edges.get(&(1, 4)).unwrap(), &[(1, 3)]);
        assert_eq!(edges.get(&(2, 3)).unwrap(), &[(1, 1)]);
        assert_eq!(edges.get(&(3, 4)).unwrap(), &[(1, 2)]);
        Ok(())
    }

    #[test]
    fn map_faces_works() -> Result<(), StrError> {
        //      8-------------11
        //     /.             /|
        //    / .            / |
        //   /  .           /  |
        //  /   .          /   |       id = 1
        // 9-------------10    |       attribute_id = 2
        // |    .         |    |
        // |    4---------|----7
        // |   /.         |   /|
        // |  / .         |  / |
        // | /  .         | /  |
        // |/   .         |/   |
        // 5--------------6    |       id = 0
        // |    .         |    |       attribute_id = 1
        // |    0---------|----3
        // |   /          |   /
        // |  /           |  /
        // | /            | /
        // |/             |/
        // 1--------------2
        let mesh = Samples::two_cubes_vertical();
        let shapes = alloc_cell_shapes(&mesh)?;
        let faces = map_faces(&mesh, &shapes);
        let mut keys: Vec<_> = faces.keys().collect();
        keys.sort();
        assert_eq!(
            keys,
            [
                &(0, 1, 2, 3),
                &(0, 1, 4, 5),
                &(0, 3, 4, 7),
                &(1, 2, 5, 6),
                &(2, 3, 6, 7),
                &(4, 5, 6, 7),
                &(4, 5, 8, 9),
                &(4, 7, 8, 11),
                &(5, 6, 9, 10),
                &(6, 7, 10, 11),
                &(8, 9, 10, 11),
            ]
        );
        assert_eq!(faces.get(&(0, 1, 2, 3)).unwrap(), &[(0, 4)]);
        assert_eq!(faces.get(&(0, 1, 4, 5)).unwrap(), &[(0, 2)]);
        assert_eq!(faces.get(&(0, 3, 4, 7)).unwrap(), &[(0, 0)]);
        assert_eq!(faces.get(&(1, 2, 5, 6)).unwrap(), &[(0, 1)]);
        assert_eq!(faces.get(&(2, 3, 6, 7)).unwrap(), &[(0, 3)]);
        assert_eq!(faces.get(&(4, 5, 6, 7)).unwrap(), &[(0, 5), (1, 4)]);
        assert_eq!(faces.get(&(4, 5, 8, 9)).unwrap(), &[(1, 2)]);
        assert_eq!(faces.get(&(4, 7, 8, 11)).unwrap(), &[(1, 0)]);
        assert_eq!(faces.get(&(5, 6, 9, 10)).unwrap(), &[(1, 1)]);
        assert_eq!(faces.get(&(6, 7, 10, 11)).unwrap(), &[(1, 3)]);
        assert_eq!(faces.get(&(8, 9, 10, 11)).unwrap(), &[(1, 5)]);

        //                       4------------7-----------10
        //                      /.           /|            |
        //                     / .          / |            |
        //                    /  .         /  |            |
        //                   /   .        /   |            |
        //                  5------------6    |            |
        //                  |    .       |`.  |            |
        //                  |    0-------|--`.3------------9
        //                  |   /        |   /`.          /
        //                  |  /         |  /   `.       /
        //                  | /          | /      `.    /
        //                  |/           |/         `. /
        //  12-----11-------1------------2------------8
        //
        let mesh = Samples::mixed_shapes_3d();
        let shapes = alloc_cell_shapes(&mesh)?;
        let faces = map_faces(&mesh, &shapes);
        let mut keys: Vec<_> = faces.keys().collect();
        keys.sort();
        assert_eq!(
            keys,
            [
                &(0, 1, 2, 3),
                &(0, 1, 4, 5),
                &(0, 3, 4, 7),
                &(1, 2, 5, 6),
                &(2, 3, 6, 7),
                &(2, 3, 6, 13),
                &(2, 3, 8, 13),
                &(2, 6, 8, 13),
                &(3, 6, 8, 13),
                &(4, 5, 6, 7),
            ]
        );
        assert_eq!(faces.get(&(0, 1, 2, 3)).unwrap(), &[(0, 4)]);
        assert_eq!(faces.get(&(0, 1, 4, 5)).unwrap(), &[(0, 2)]);
        assert_eq!(faces.get(&(0, 3, 4, 7)).unwrap(), &[(0, 0)]);
        assert_eq!(faces.get(&(1, 2, 5, 6)).unwrap(), &[(0, 1)]);
        assert_eq!(faces.get(&(2, 3, 6, 7)).unwrap(), &[(0, 3)]);
        assert_eq!(faces.get(&(2, 3, 6, 13)).unwrap(), &[(1, 0)]);
        assert_eq!(faces.get(&(2, 3, 8, 13)).unwrap(), &[(1, 2)]);
        assert_eq!(faces.get(&(2, 6, 8, 13)).unwrap(), &[(1, 1)]);
        assert_eq!(faces.get(&(3, 6, 8, 13)).unwrap(), &[(1, 3)]);
        assert_eq!(faces.get(&(4, 5, 6, 7)).unwrap(), &[(0, 5)]);
        Ok(())
    }
}

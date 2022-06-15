use crate::mesh::{EdgeKey, MapEdge2dToCells, Mesh};
use russell_lab::sort2;
use std::collections::HashMap;

/// Extracts all 2d edges (internal and boundary)
///
/// # Output
///
/// * Returns a map relating edge keys to `Vec<(cell_id, e)>` where:
///     - `cell_id` -- is the id of the cell sharing the edge
///     - `e` -- is the cell's local edge index
///
/// # Panics
///
/// 1. It panics if `shapes.len() != mesh.cells.len()` (i.e., the shapes vector and the mesh must be compatible)
/// 2. It panics if `mesh.ndim != 2` (i.e., this function works in 2D only)
pub(crate) fn extract_all_2d_edges(mesh: &Mesh) -> MapEdge2dToCells {
    assert_eq!(mesh.ndim, 2);
    let mut edges = HashMap::new();
    mesh.cells.iter().for_each(|cell| {
        if cell.kind.ndim() == 2 {
            for e in 0..cell.kind.nedge() {
                let mut edge_key: EdgeKey = (
                    cell.points[cell.kind.edge_node_id(e, 0)],
                    cell.points[cell.kind.edge_node_id(e, 1)],
                );
                sort2(&mut edge_key);
                let data = edges.entry(edge_key).or_insert(Vec::new());
                data.push((cell.id, e));
            }
        }
    });
    edges
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::extract_all_2d_edges;
    use crate::mesh::Samples;

    #[test]
    fn extract_all_2d_edges_works() {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let edges = extract_all_2d_edges(&mesh);
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
    }

    #[test]
    fn extract_all_2d_edges_mixed_works() {
        //           4---------3
        //           |         |
        //           |   [1]   |
        //           |         |
        //  0--------1---------2
        let mesh = Samples::mixed_shapes_2d();
        let edges = extract_all_2d_edges(&mesh);
        let mut keys: Vec<_> = edges.keys().collect();
        keys.sort();
        assert_eq!(keys, [&(1, 2), &(1, 4), &(2, 3), &(3, 4)]);
        assert_eq!(edges.get(&(1, 2)).unwrap(), &[(1, 0)]);
        assert_eq!(edges.get(&(1, 4)).unwrap(), &[(1, 3)]);
        assert_eq!(edges.get(&(2, 3)).unwrap(), &[(1, 1)]);
        assert_eq!(edges.get(&(3, 4)).unwrap(), &[(1, 2)]);
    }
}

use crate::mesh::{EdgeKey, EdgesCellsMap2D, Mesh};
use crate::shapes::Shape;
use crate::StrError;
use russell_lab::sort2;
use std::collections::HashMap;

/// Extracts all edges (internal and boundary) in 2D
///
/// # Input
///
/// * `mesh` -- the Mesh
/// * `shapes` -- the shapes of cells (len == cells.len())
///
/// # Output
///
/// * Returns a map relating edge keys to `Vec<(cell_id, e)>` where:
///     - `cell_id` -- the id of the cell sharing the edge
///     - `e` -- is the cell's local edge index
#[inline]
pub(crate) fn extract_all_edges_2d(mesh: &Mesh, shapes: &Vec<Shape>) -> Result<EdgesCellsMap2D, StrError> {
    if mesh.space_ndim != 2 {
        return Err("this function works in 2D only");
    }
    let mut edges: EdgesCellsMap2D = HashMap::new();
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
    Ok(edges)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::extract_all_edges_2d;
    use crate::mesh::{allocate_shapes, Mesh, Samples};
    use crate::StrError;

    #[test]
    fn capture_some_wrong_input() {
        let mesh = Mesh {
            space_ndim: 1,
            points: Vec::new(),
            cells: Vec::new(),
        };
        assert_eq!(
            extract_all_edges_2d(&mesh, &Vec::new()).err(),
            Some("this function works in 2D only")
        );
    }

    #[test]
    fn extract_all_edges_2d_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let shapes = allocate_shapes(&mesh)?;
        let edges = extract_all_edges_2d(&mesh, &shapes)?;
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
        Ok(())
    }

    #[test]
    fn extract_all_edges_2d_mixed_works() -> Result<(), StrError> {
        //           4---------3
        //           |         |
        //           |   [1]   |
        //           |         |
        //  0--------1---------2
        let mesh = Samples::mixed_shapes_2d();
        let shapes = allocate_shapes(&mesh)?;
        let edges = extract_all_edges_2d(&mesh, &shapes)?;
        let mut keys: Vec<_> = edges.keys().collect();
        keys.sort();
        assert_eq!(keys, [&(1, 2), &(1, 4), &(2, 3), &(3, 4)]);
        assert_eq!(edges.get(&(1, 2)).unwrap(), &[(1, 0)]);
        assert_eq!(edges.get(&(1, 4)).unwrap(), &[(1, 3)]);
        assert_eq!(edges.get(&(2, 3)).unwrap(), &[(1, 1)]);
        assert_eq!(edges.get(&(3, 4)).unwrap(), &[(1, 2)]);
        Ok(())
    }
}

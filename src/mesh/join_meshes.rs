use super::{Cell, Features, Mesh, Point};
use crate::mesh::Extract;
use crate::util::GridSearch;
use crate::StrError;

/// Joins two meshes by comparing the coordinates on the boundary of the first mesh
///
/// **Note:** The meshes must have the same space_ndim.
///
/// **Important:** This function does not guarantee the "mesh compatibility" requirements
/// for finite element analyses.
pub fn join_meshes(a: &Mesh, b: &Mesh) -> Result<Mesh, StrError> {
    // check
    if a.ndim != b.ndim {
        return Err("meshes must have the same ndim");
    }

    // find the boundary of mesh A
    let (_, _, boundary_a) = Features::new(a, Extract::Boundary);

    // allocate and prepare a GridSearch for mesh A
    let mut grid_a = GridSearch::new(&boundary_a.min, &boundary_a.max, None, None, None)?;
    for m in 0..a.points.len() {
        grid_a.insert(a.points[m].id, &a.points[m].coords)?;
    }

    // create a new mesh with all mesh A data
    let mut mesh = a.clone();

    // renumber the points of mesh B and add them to the new mesh (if not present yet)
    let mut new_point_id = a.points.len();
    let mut map_old_to_new_point_id_b = vec![0; b.points.len()];
    for m in 0..b.points.len() {
        let x = &b.points[m].coords;
        let maybe_point_id_a = if grid_a.is_outside(x) { None } else { grid_a.find(x)? };
        let id = match maybe_point_id_a {
            Some(point_id_a) => point_id_a,
            None => {
                mesh.points.push(Point {
                    id: new_point_id,
                    coords: x.clone(),
                });
                new_point_id += 1;
                new_point_id - 1
            }
        };
        map_old_to_new_point_id_b[m] = id;
    }

    // insert cells of mesh B into new mesh
    let mut new_cell_id = a.cells.len();
    for cell in &b.cells {
        mesh.cells.push(Cell {
            id: new_cell_id,
            attribute_id: cell.attribute_id,
            kind: cell.kind,
            points: cell.points.iter().map(|id| map_old_to_new_point_id_b[*id]).collect(),
        });
        new_cell_id += 1;
    }
    Ok(mesh)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::join_meshes;
    use crate::mesh::{check_ids_and_kind, check_jacobian, Samples};
    use crate::StrError;

    #[test]
    fn join_meshes_handles_errors() {
        let a = Samples::two_qua4();
        let b = Samples::two_hex8();
        assert_eq!(join_meshes(&a, &b).err(), Some("meshes must have the same ndim"));
    }

    #[test]
    fn join_meshes_works_2d() -> Result<(), StrError> {
        //          [#] indicates id
        //      y   (#) indicates attribute_id
        //      ↑
        // 1.0  3-----------2-----------5
        //      |           |           |
        //      |    [0]    |    [1]    |
        //      |    (1)    |    (2)    |
        //      |           |           |
        // 0.0  0-----------1-----------4  → x
        //
        //                  +
        //
        // 1.0  3-----------2-----------5
        //      |           |           |
        //      |    [0]    |    [1]    |
        //      |    (1)    |    (2)    |
        //      |           |           |
        // 0.0  0-----------1-----------4  → x
        //     0.0         1.0         2.0
        let a = Samples::two_qua4();
        let mut b = Samples::two_qua4();

        // shift B-mesh up
        for m in 0..b.points.len() {
            b.points[m].coords[1] += 1.0;
        }

        //      y
        //      ↑
        // 2.0  7-----------6-----------8
        //      |           |           |
        //      |    [2]    |    [3]    |
        //      |    (1)    |    (2)    |
        //      |           |           |
        // 1.0  3-----------2-----------5
        //      |           |           |
        //      |    [0]    |    [1]    |
        //      |    (1)    |    (2)    |
        //      |           |           |
        // 0.0  0-----------1-----------4  → x
        //     0.0         1.0         2.0
        let mesh = join_meshes(&a, &b)?;
        check_ids_and_kind(&mesh)?;
        check_jacobian(&mesh)?;
        assert_eq!(mesh.points[0].coords, &[0.0, 0.0]);
        assert_eq!(mesh.points[1].coords, &[1.0, 0.0]);
        assert_eq!(mesh.points[2].coords, &[1.0, 1.0]);
        assert_eq!(mesh.points[3].coords, &[0.0, 1.0]);
        assert_eq!(mesh.points[4].coords, &[2.0, 0.0]);
        assert_eq!(mesh.points[5].coords, &[2.0, 1.0]);
        assert_eq!(mesh.points[6].coords, &[1.0, 2.0]);
        assert_eq!(mesh.points[7].coords, &[0.0, 2.0]);
        assert_eq!(mesh.points[8].coords, &[2.0, 2.0]);
        assert_eq!(mesh.cells[0].points, &[0, 1, 2, 3]);
        assert_eq!(mesh.cells[1].points, &[1, 4, 5, 2]);
        assert_eq!(mesh.cells[2].points, &[3, 2, 6, 7]);
        assert_eq!(mesh.cells[3].points, &[2, 5, 8, 6]);
        assert_eq!(mesh.cells[0].attribute_id, 1);
        assert_eq!(mesh.cells[1].attribute_id, 2);
        assert_eq!(mesh.cells[2].attribute_id, 1);
        assert_eq!(mesh.cells[3].attribute_id, 2);
        Ok(())
    }

    #[test]
    fn join_meshes_works_3d() -> Result<(), StrError> {
        //       8-------------11  2.0          8-------------11  2.0
        //      /.             /|              /.             /|
        //     / .            / |             / .            / |
        //    /  .           /  |            /  .           /  |
        //   /   .          /   |           /   .          /   |
        //  9-------------10    |          9-------------10    |
        //  |    .         |    |          |    .         |    |
        //  |    4---------|----7  1.0     |    4---------|----7  1.0
        //  |   /. [1]     |   /|          |   /. [1]     |   /|
        //  |  / . (2)     |  / |      +   |  / . (2)     |  / |
        //  | /  .         | /  |          | /  .         | /  |
        //  |/   .         |/   |          |/   .         |/   |
        //  5--------------6    |          5--------------6    |          z
        //  |    .         |    |          |    .         |    |          ↑
        //  |    0---------|----3  0.0     |    0---------|----3  0.0     o → y
        //  |   /  [0]     |   /           |   /  [0]     |   /          ↙
        //  |  /   (1)     |  /            |  /   (1)     |  /          x
        //  | /            | /             | /            | /
        //  |/             |/              |/             |/
        //  1--------------2   1.0         1--------------2   1.0
        let a = Samples::two_hex8();
        let mut b = Samples::two_hex8();

        // shift mesh B along y (to the right)
        for m in 0..b.points.len() {
            b.points[m].coords[1] += 1.0;
        }

        //       8-------------11-------------17  2.0
        //      /.             /.             /|
        //     / .            / .            / |
        //    /  .           /  .           /  |
        //   /   .          /   .          /   |
        //  9=============10=============16    |
        //  |    .         |    .         |    |
        //  |    4---------|----7---------|---15  1.0
        //  |   /. [1]     |   /. [3]     |   /|
        //  |  / . (2)     |  / . (2)     |  / |
        //  | /  .         | /  .         | /  |
        //  |/   .         |/   .         |/   |
        //  5--------------6-------------14    |          z
        //  |    .         |    .         |    |          ↑
        //  |    0---------|----3---------|---13  0.0     o → y
        //  |   /  [0]     |   /  [2]     |   /          ↙
        //  |  /   (1)     |  /   (1)     |  /          x
        //  | /            | /            | /
        //  |/             |/             |/
        //  1--------------2-------------12   1.0
        // 0.0            1.0            2.0
        let mesh = join_meshes(&a, &b)?;
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 18);
        assert_eq!(mesh.cells.len(), 4);
        check_ids_and_kind(&mesh)?;
        check_jacobian(&mesh)?;

        let sample = Samples::four_hex8();
        for i in 0..mesh.cells.len() {
            assert_eq!(mesh.cells[i].id, sample.cells[i].id);
            assert_eq!(mesh.cells[i].points, sample.cells[i].points);
        }
        for i in 0..mesh.points.len() {
            assert_eq!(mesh.points[i].id, sample.points[i].id);
            assert_eq!(mesh.points[i].coords, sample.points[i].coords);
        }
        Ok(())
    }
}

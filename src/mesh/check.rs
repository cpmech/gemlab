use super::{allocate_state, Mesh};
use crate::shapes::Shape;
use crate::StrError;

/// Checks if the IDs of points and cells are sequential
///
/// This function checks that:
///
/// * the index of a point in the points vector matches the id of the point
/// * the index of a cell in the cells vector matches the id of the cell
pub fn check_ids(mesh: &Mesh) -> Result<(), StrError> {
    for i in 0..mesh.points.len() {
        if mesh.points[i].id != i {
            return Err("invalid point id found; ids must be sequential");
        }
    }
    for i in 0..mesh.cells.len() {
        if mesh.cells[i].id != i {
            return Err("invalid cell id found; ids must be sequential");
        }
    }
    Ok(())
}

/// Checks if the determinant of the Jacobian of all cells are non-negative
pub fn check_jacobian(mesh: &Mesh) -> Result<(), StrError> {
    let ksi = [0.0, 0.0, 0.0];
    for cell in &mesh.cells {
        let shape = Shape::new(mesh.space_ndim, cell.geo_ndim, cell.points.len())?;
        let mut state = allocate_state(mesh, cell.geo_ndim, &cell.points)?;
        let det_jac = shape.calc_jacobian(&mut state, &ksi)?;
        if det_jac < 0.0 {
            return Err("negative determinant of Jacobian found");
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{check_ids, check_jacobian};
    use crate::mesh::{Cell, Mesh, Point};

    #[test]
    fn check_ids_works() {
        //  3-----------2-----------5
        //  |           |           |
        //  |    [0]    |    [1]    |
        //  |    (1)    |    (2)    |
        //  |           |           |
        //  0-----------1-----------4
        #[rustfmt::skip]
        let mut mesh = Mesh {
            space_ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, geo_ndim: 2, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute_id: 2, geo_ndim: 2, points: vec![1, 4, 5, 2] },
            ],
        };
        assert_eq!(check_ids(&mesh).expect("should not fail"), ());

        mesh.points[0].id = 6;
        assert_eq!(
            check_ids(&mesh).err(),
            Some("invalid point id found; ids must be sequential")
        );
        mesh.points[0].id = 0;

        assert_eq!(check_ids(&mesh).expect("should not fail"), ());

        mesh.cells[0].id = 2;
        assert_eq!(
            check_ids(&mesh).err(),
            Some("invalid cell id found; ids must be sequential")
        );
        mesh.cells[0].id = 0;
    }

    #[test]
    fn check_jacobian_works() {
        //  3-----------2-----------5
        //  |           |           |
        //  |    [0]    |    [1]    |
        //  |    (1)    |    (2)    |
        //  |           |           |
        //  0-----------1-----------4
        #[rustfmt::skip]
        let mut mesh = Mesh {
            space_ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, geo_ndim: 2, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute_id: 2, geo_ndim: 2, points: vec![1, 4, 5, 2] },
            ],
        };
        assert_eq!(check_jacobian(&mesh).expect("should not fail"), ());

        mesh.cells[0].points[2] = 3;
        mesh.cells[0].points[3] = 2;
        assert_eq!(
            check_jacobian(&mesh).err(),
            Some("cannot compute inverse due to zero determinant")
        );

        mesh.cells[0].points = vec![0, 3, 2, 1];
        assert_eq!(
            check_jacobian(&mesh).err(),
            Some("negative determinant of Jacobian found")
        );
    }
}

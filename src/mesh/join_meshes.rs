use super::{allocate_cell_shapes, Cell, Features, Mesh, Point};
use crate::mesh::Extract;
use crate::util::{GridSearch, GsNdiv, GsTol};
use crate::StrError;

/// Joins two meshes by comparing the coordinates on the boundary of the first mesh
///
/// **Note:** The meshes must have the same space_ndim.
///
/// **Important:** This function does not guarantee the "mesh compatibility" requirements
/// for finite element analyses.
pub fn join_meshes(a: &Mesh, b: &Mesh) -> Result<Mesh, StrError> {
    // check
    if a.space_ndim != b.space_ndim {
        return Err("meshes must have the same space_ndim");
    }

    // find the boundary of mesh A
    let shapes_a = allocate_cell_shapes(&a)?;
    let (_, _, boundary_a) = Features::new(a, &shapes_a, Extract::Boundary);

    // allocate and prepare a GridSearch for mesh A
    let mut min_a = boundary_a.min.clone();
    let mut max_a = boundary_a.max.clone();
    const PCT: f64 = 1.0 / 100.0;
    for i in 0..a.space_ndim {
        let del = max_a[i] - min_a[i];
        min_a[i] -= PCT * del;
        max_a[i] += PCT * del;
    }
    let mut grid_a = GridSearch::new(&min_a, &max_a, GsNdiv::Default, GsTol::Default)?;
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
        let mut in_grid_a = true;
        for i in 0..a.space_ndim {
            if x[i] < min_a[i] || x[i] > max_a[i] {
                in_grid_a = false;
                break;
            }
        }
        let maybe_point_id_a = if in_grid_a { grid_a.find(x)? } else { None };
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
            geo_ndim: cell.geo_ndim,
            points: cell.points.iter().map(|id| map_old_to_new_point_id_b[*id]).collect(),
        });
        new_cell_id += 1;
    }
    Ok(mesh)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}

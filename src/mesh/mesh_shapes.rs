use super::{EdgeKey, FaceKey, Mesh};
use crate::{shapes::Shape, StrError};
use std::collections::HashMap;

pub struct MeshShapes {
    pub cells: Vec<Shape>,
    pub boundary_edges: HashMap<EdgeKey, Shape>,
    pub boundary_faces: HashMap<FaceKey, Shape>,
}

impl MeshShapes {
    pub fn new(mesh: &Mesh) -> Result<Self, StrError> {
        let mut fem_mesh = MeshShapes {
            cells: Vec::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
        };
        for cell in &mesh.cells {
            let nnode = cell.points.len();
            let mut shape = Shape::new(mesh.space_ndim, cell.geo_ndim, nnode)?;
            for m in 0..nnode {
                let point_id = cell.points[m];
                for j in 0..mesh.space_ndim {
                    shape.set_node(m, j, mesh.points[point_id].coords[j])?;
                }
            }
            fem_mesh.cells.push(shape);
        }
        for (edge_key, edge) in &mesh.boundary_edges {
            let nnode = edge.points.len();
            let mut shape = Shape::new(mesh.space_ndim, 1, nnode)?;
            for m in 0..nnode {
                let point_id = edge.points[m];
                for j in 0..mesh.space_ndim {
                    shape.set_node(m, j, mesh.points[point_id].coords[j])?;
                }
            }
            fem_mesh.boundary_edges.insert(*edge_key, shape);
        }
        for (face_key, face) in &mesh.boundary_faces {
            let nnode = face.points.len();
            let mut shape = Shape::new(mesh.space_ndim, 1, nnode)?;
            for m in 0..nnode {
                let point_id = face.points[m];
                for j in 0..mesh.space_ndim {
                    shape.set_node(m, j, mesh.points[point_id].coords[j])?;
                }
            }
            fem_mesh.boundary_faces.insert(*face_key, shape);
        }
        Ok(fem_mesh)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{mesh::Mesh, StrError};

    use super::MeshShapes;

    #[test]
    fn from_text_and_strings_work() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mesh = Mesh::from_text(
            r"# header
            # ndim npoint ncell
                 2      6     2
            
            # points
            # id   x   y
               0 0.0 0.0
               1 1.0 0.0
               2 1.0 1.0
               3 0.0 1.0
               4 2.0 0.0
               5 2.0 1.0
            
            # cells
            # id att geo_ndim npoint point_ids...
               0   1        2      4 0 1 2 3
               1   0        2      4 1 4 5 2",
        )?;

        let mesh_shapes = MeshShapes::new(&mesh)?;
        assert_eq!(
            format!("{}", mesh_shapes.cells[0].coords_transp),
            "┌         ┐\n\
             │ 0 1 1 0 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        assert_eq!(
            format!("{}", mesh_shapes.cells[1].coords_transp),
            "┌         ┐\n\
             │ 1 2 2 1 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        assert_eq!(
            format!("{}", mesh_shapes.boundary_edges.get(&(0, 1)).unwrap().coords_transp),
            "┌     ┐\n\
             │ 0 1 │\n\
             │ 0 0 │\n\
             └     ┘"
        );
        assert_eq!(
            format!("{}", mesh_shapes.boundary_edges.get(&(1, 4)).unwrap().coords_transp),
            "┌     ┐\n\
             │ 1 2 │\n\
             │ 0 0 │\n\
             └     ┘"
        );
        assert_eq!(
            format!("{}", mesh_shapes.boundary_edges.get(&(2, 3)).unwrap().coords_transp),
            "┌     ┐\n\
             │ 1 0 │\n\
             │ 1 1 │\n\
             └     ┘"
        );
        assert_eq!(mesh_shapes.boundary_faces.len(), 0);

        Ok(())
    }
}

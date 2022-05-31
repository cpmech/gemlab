use crate::mesh::{Edge, EdgeKey, Extract, Face, FacesCellsMap3D, Features, Mesh};
use crate::shapes::Shape;
use russell_lab::sort2;
use std::collections::{HashMap, HashSet};

/// Extracts mesh features in 3D
pub(crate) fn extract_features_3d(
    mesh: &Mesh,
    shapes: &Vec<Shape>,
    faces: &FacesCellsMap3D,
    extract: Extract,
) -> Features {
    assert_eq!(mesh.space_ndim, 3);

    // output
    let mut features = Features {
        points: HashSet::new(),
        edges: HashMap::new(),
        faces: HashMap::new(),
        min: vec![f64::MAX; mesh.space_ndim],
        max: vec![f64::MIN; mesh.space_ndim],
    };

    // sort face keys just so the next loop is deterministic
    let mut face_keys: Vec<_> = faces.keys().collect();
    face_keys.sort();

    // loop over all faces
    for face_key in face_keys {
        let shared_by = faces.get(face_key).unwrap();

        // accept feature?
        let accept = match extract {
            Extract::All => true,
            Extract::Boundary => shared_by.len() == 1, // boundary face; with only one shared cell
            Extract::Interior => shared_by.len() != 1, // interior face; shared by multiple cells
        };
        if !accept {
            continue;
        }

        // cell and face
        let (cell_id, f) = shared_by[0];
        let cell = &mesh.cells[cell_id];
        let shape = &shapes[cell_id];
        let mut face = Face {
            points: vec![0; shape.face_nnode],
        };

        // process points on face
        for i in 0..shape.face_nnode {
            face.points[i] = cell.points[shape.face_node_id(f, i)];
            features.points.insert(face.points[i]);
            for j in 0..mesh.space_ndim {
                features.min[j] = f64::min(features.min[j], mesh.points[face.points[i]].coords[j]);
                features.max[j] = f64::max(features.max[j], mesh.points[face.points[i]].coords[j]);
            }
        }

        // loop over all edges on face
        let face_shape = Shape::new(mesh.space_ndim, 2, shape.face_nnode).unwrap(); // should not fail
        for e in 0..face_shape.nedge {
            // define edge key (sorted point ids)
            let mut edge_key: EdgeKey = (
                face.points[face_shape.edge_node_id(e, 0)],
                face.points[face_shape.edge_node_id(e, 1)],
            );
            sort2(&mut edge_key);

            // skip already handled edge
            if features.edges.contains_key(&edge_key) {
                continue;
            }

            // new edge
            let mut edge = Edge {
                points: vec![0; face_shape.edge_nnode],
            };
            for i in 0..face_shape.edge_nnode {
                edge.points[i] = face.points[face_shape.edge_node_id(e, i)];
            }
            features.edges.insert(edge_key, edge);
        }

        // new face
        features.faces.insert(*face_key, face);
    }
    features
}

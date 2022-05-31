use crate::mesh::{Edge, Extract, Features, MapEdge2dToCells, Mesh};
use crate::shapes::Shape;
use std::collections::{HashMap, HashSet};

/// Extracts mesh features in 2D
pub(crate) fn extract_features_2d(
    mesh: &Mesh,
    shapes: &Vec<Shape>,
    edges: &MapEdge2dToCells,
    extract: Extract,
) -> Features {
    assert_eq!(mesh.space_ndim, 2);

    // output
    let mut features = Features {
        points: HashSet::new(),
        edges: HashMap::new(),
        faces: HashMap::new(),
        min: vec![f64::MAX; mesh.space_ndim],
        max: vec![f64::MIN; mesh.space_ndim],
    };

    // loop over all edges
    for (edge_key, shared_by) in edges {
        // accept feature?
        let accept = match extract {
            Extract::All => true,
            Extract::Boundary => shared_by.len() == 1, // boundary edge; with only one shared cell
            Extract::Interior => shared_by.len() != 1, // interior edge; shared by multiple cells
        };
        if !accept {
            continue;
        }

        // cell and edge
        let (cell_id, e) = shared_by[0];
        let cell = &mesh.cells[cell_id];
        let shape = &shapes[cell_id];
        let mut edge = Edge {
            points: vec![0; shape.edge_nnode],
        };

        // process points on edge
        for i in 0..shape.edge_nnode {
            edge.points[i] = cell.points[shape.edge_node_id(e, i)];
            features.points.insert(edge.points[i]);
            for j in 0..mesh.space_ndim {
                features.min[j] = f64::min(features.min[j], mesh.points[edge.points[i]].coords[j]);
                features.max[j] = f64::max(features.max[j], mesh.points[edge.points[i]].coords[j]);
            }
        }

        // new edge
        features.edges.insert(*edge_key, edge);
    }
    features
}

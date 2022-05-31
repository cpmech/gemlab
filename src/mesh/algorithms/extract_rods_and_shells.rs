use crate::mesh::{Features, Mesh};

/// Extracts points of rods (2D or 3D) or shells (3D)
#[inline]
pub(crate) fn extract_rods_and_shells(mesh: &Mesh, features: &mut Features) {
    mesh.cells.iter().for_each(|cell| {
        if cell.geo_ndim == 1 || (cell.geo_ndim == 2 && mesh.space_ndim == 3) {
            cell.points.iter().for_each(|id| {
                features.points.insert(*id);
                for j in 0..mesh.space_ndim {
                    features.min[j] = f64::min(features.min[j], mesh.points[*id].coords[j]);
                    features.max[j] = f64::max(features.max[j], mesh.points[*id].coords[j]);
                }
            });
        }
    });
}

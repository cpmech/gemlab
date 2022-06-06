use crate::mesh::{Features, Mesh};

/// Extracts points of rods (2D or 3D) or shells (3D)
#[inline]
pub(crate) fn extract_rods_and_shells(mesh: &Mesh, features: &mut Features) {
    mesh.cells.iter().for_each(|cell| {
        let geo_ndim = cell.kind.ndim();
        if geo_ndim == 1 || (geo_ndim == 2 && mesh.ndim == 3) {
            cell.points.iter().for_each(|id| {
                features.points.insert(*id);
                for j in 0..mesh.ndim {
                    features.min[j] = f64::min(features.min[j], mesh.points[*id].coords[j]);
                    features.max[j] = f64::max(features.max[j], mesh.points[*id].coords[j]);
                }
            });
        }
    });
}

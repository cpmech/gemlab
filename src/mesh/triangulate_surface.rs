use super::AsCell;
use crate::mesh::{GeoClass, Mesh, PointId};
use crate::shapes::{GeoKind, Scratchpad};
use russell_lab::Vector;
use std::collections::HashMap;

/// Holds the result of a triangulation
pub struct Triangulation {
    /// X coordinates of the points
    pub xx: Vec<f64>,

    /// Y coordinates of the points
    pub yy: Vec<f64>,

    /// Z coordinates of the points
    pub zz: Vec<f64>,

    /// Triangles (connectivity) defined by point indices
    pub triangles: Vec<Vec<PointId>>,
}

impl Triangulation {
    /// Triangulates a surface defined by Cells or Faces (AsCell) in 2D or 3D
    ///
    /// **Note:** Only Tri and Qua cells will be triangulated; other cell classes will be ignored.
    /// Therefore, the result may be empty if no Tri or Qua cells are found on the surface.
    pub fn from_surface<T: AsCell>(mesh: &Mesh, surface: &[&T]) -> Self {
        // constants
        let ndim = mesh.ndim;

        // auxiliary
        let mut old_point_id_to_new_point_id: HashMap<PointId, PointId> = HashMap::new();
        let mut pads: HashMap<GeoKind, Scratchpad> = HashMap::new();
        let mut x_work = Vector::new(ndim);
        let mut connectivity = vec![0; 3]; // 3 nodes per triangle

        // results
        let mut res = Triangulation {
            xx: Vec::new(),
            yy: Vec::new(),
            zz: Vec::new(),
            triangles: Vec::new(),
        };

        // loop over all surface cells
        for cell in surface {
            // skip non-Tri/Qua cells
            let kind = cell.kind();
            let class = kind.class();
            let tri_or_qua = class == GeoClass::Tri || class == GeoClass::Qua;
            if !tri_or_qua {
                continue;
            }

            // constants
            let points = cell.points();
            let nnode = kind.nnode();
            let extra_nnode = kind.triangulate_extra_nnode();
            let mut extra_points = vec![usize::MAX; extra_nnode];

            // set scratchpad
            let pad = pads.entry(kind).or_insert_with(|| Scratchpad::new(ndim, kind).unwrap());
            mesh.set_pad(pad, points);

            // loop over triangles and add them to the list of triangles
            pad.triangulate(&mut x_work, |_, i, m, x| {
                let p = if m < nnode {
                    // existing point
                    let p_old = points[m];
                    old_point_id_to_new_point_id
                        .entry(p_old)
                        .or_insert_with(|| {
                            let p_new = res.xx.len();
                            res.xx.push(x[0]);
                            res.yy.push(x[1]);
                            if ndim == 3 {
                                res.zz.push(x[2]);
                            }
                            p_new
                        })
                        .clone()
                } else {
                    // extra point
                    let k = m - nnode;
                    if extra_points[k] == usize::MAX {
                        let p_new = res.xx.len();
                        extra_points[k] = p_new;
                        res.xx.push(x[0]);
                        res.yy.push(x[1]);
                        if ndim == 3 {
                            res.zz.push(x[2]);
                        }
                        p_new
                    } else {
                        extra_points[k]
                    }
                };
                connectivity[i] = p;
                if i == 2 {
                    // last corner
                    res.triangles.push(connectivity.clone());
                }
            })
            .unwrap();
        }

        // done
        res
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Triangulation;
    use crate::mesh::Samples;

    #[test]
    fn triangulate_surface_works_1() {
        let mesh = Samples::one_qua4();
        let surface: Vec<_> = mesh.cells.iter().collect();
        let res = Triangulation::from_surface(&mesh, &surface);
        assert_eq!(res.xx.len(), 4);
        assert_eq!(res.yy.len(), 4);
        assert_eq!(res.zz.len(), 0);
        assert_eq!(res.triangles.len(), 2);
    }
}

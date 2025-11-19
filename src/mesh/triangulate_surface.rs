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
    use plotpy::{Canvas, Plot};

    const SAVE_FIGURE: bool = false;

    #[test]
    fn triangulate_surface_works_1() {
        let mesh = Samples::one_qua4();
        let surface: Vec<_> = mesh.cells.iter().collect();
        let res = Triangulation::from_surface(&mesh, &surface);
        if SAVE_FIGURE {
            let mut canvas = Canvas::new();
            canvas.draw_triangles(&res.xx, &res.yy, &res.triangles);
            let mut plot = Plot::new();
            plot.add(&canvas)
                .set_equal_axes(true)
                .save("/tmp/gemlab/triangulate_surface_works_1.svg")
                .unwrap();
        }
        assert_eq!(res.xx.len(), 4);
        assert_eq!(res.yy.len(), 4);
        assert_eq!(res.zz.len(), 0);
        assert_eq!(res.triangles.len(), 2);
    }

    #[test]
    fn triangulate_surface_works_2() {
        let mesh = Samples::qua8_tri6_lin2();
        let surface: Vec<_> = mesh.cells.iter().collect();
        let res = Triangulation::from_surface(&mesh, &surface);
        if SAVE_FIGURE {
            let mut canvas = Canvas::new();
            canvas.draw_triangles(&res.xx, &res.yy, &res.triangles);
            let mut plot = Plot::new();
            plot.add(&canvas)
                .set_equal_axes(true)
                .save("/tmp/gemlab/triangulate_surface_works_2.svg")
                .unwrap();
        }
        let npoint = mesh.points.len() + 1; // extra point at the center of Qua8
        assert_eq!(res.xx.len(), npoint);
        assert_eq!(res.yy.len(), npoint);
        assert_eq!(res.zz.len(), 0);
        assert_eq!(res.triangles.len(), 8 + 4); // 8 from Qua8 and 4 from Tri6
    }

    #[test]
    fn triangulate_surface_works_3() {
        let mesh = Samples::qua8_tri6_lin2_three_dimensional();
        let surface: Vec<_> = mesh.cells.iter().collect();
        let res = Triangulation::from_surface(&mesh, &surface);
        if SAVE_FIGURE {
            let mut canvas = Canvas::new();
            canvas.draw_triangles_3d(&res.xx, &res.yy, &res.zz, &res.triangles);
            let mut plot = Plot::new();
            plot.add(&canvas)
                .set_equal_axes(true)
                .set_camera(30.0, 90.0)
                .set_figure_size_points(600.0, 600.0)
                .save("/tmp/gemlab/triangulate_surface_works_3.svg")
                .unwrap();
        }
        let npoint = mesh.points.len() + 1; // extra point at the center of Qua8
        assert_eq!(res.xx.len(), npoint);
        assert_eq!(res.yy.len(), npoint);
        assert_eq!(res.zz.len(), npoint);
        assert_eq!(res.triangles.len(), 8 + 4); // 8 from Qua8 and 4 from Tri6
    }

    #[test]
    fn triangulate_surface_works_4() {
        let mesh = Samples::block_2d_four_qua12();
        let surface: Vec<_> = mesh.cells.iter().collect();
        let res = Triangulation::from_surface(&mesh, &surface);
        if SAVE_FIGURE {
            let mut canvas = Canvas::new();
            canvas.draw_triangles(&res.xx, &res.yy, &res.triangles);
            let mut plot = Plot::new();
            plot.add(&canvas)
                .set_equal_axes(true)
                .save("/tmp/gemlab/triangulate_surface_works_4.svg")
                .unwrap();
        }
        let ncell = 4;
        let ntriangle = ncell * 18; // 4 Qua12, each triangulated into 2*9 triangles
        let npoint = mesh.points.len() + 4 * ncell; // 4 extra points per Qua12
        assert_eq!(res.xx.len(), npoint);
        assert_eq!(res.yy.len(), npoint);
        assert_eq!(res.zz.len(), 0);
        assert_eq!(res.triangles.len(), ntriangle);
    }
}

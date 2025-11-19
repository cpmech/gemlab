use super::{GeoClass, Scratchpad};
use crate::StrError;
use russell_lab::Vector;

impl Scratchpad {
    /// Triangulates a Tri or Qua shape
    ///
    /// This function connects the existing nodes of a GeoKind to make triangles.
    /// For some shapes, a few extra nodes are defined.
    ///
    /// The results are available via the callback function: `f(t, i, m, x)` where:
    ///
    /// * `t` -- is the index of the triangle in `0..kind.triangulate_ntriangle()`
    /// * `i` -- is the corner index in `0..3`
    /// * `m` -- is the node number in `0..(kind.nnode()+kind.triangulate_extra_nnode())`
    /// * `x` -- holds the real coordinates of the node
    pub fn triangulate<F>(&mut self, mut f: F) -> Result<(), StrError>
    where
        F: FnMut(usize, usize, usize, &Vector),
    {
        let kind = self.kind;
        let class = kind.class();
        let tri_or_qua = class == GeoClass::Tri || class == GeoClass::Qua;
        if !tri_or_qua {
            return Err("triangulate works with Tri and Qua classes only");
        }
        let space_ndim = self.get_space_ndim();
        let mut x = Vector::new(space_ndim);
        for t in 0..kind.triangulate_ntriangle() {
            for i in 0..3 {
                let m = kind.triangulate_triangle_nodes(t, i);
                let ksi = if m >= kind.nnode() {
                    let k = m - kind.nnode();
                    kind.triangulate_extra_coords(k)
                } else {
                    kind.reference_coords(m)
                };
                self.calc_coords(&mut x, ksi)?;
                f(t, i, m, &x);
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::geometry::{Point2d, Triangle2d};
    use crate::shapes::scratchpad_testing::aux;
    use crate::shapes::{GeoClass, GeoKind};
    use plotpy::{Canvas, Plot};
    use russell_lab::approx_eq;
    use russell_lab::math::PI;

    const SAVE_FIGURE: bool = false;

    #[test]
    fn pad_triangulate_works() {
        let ninety = PI / 2.0;
        let forty_five = PI / 4.0;
        for kind in GeoKind::VALUES {
            let ntriangle = kind.triangulate_ntriangle();
            if ntriangle > 0 {
                let edge_nnode = kind.edge_nnode();
                let n_quads = (edge_nnode - 1) * (edge_nnode - 1);
                if kind.class() == GeoClass::Tri {
                    assert_eq!(ntriangle, n_quads);
                } else {
                    assert_eq!(ntriangle, n_quads * 2);
                }
                let expected_area = if kind.class() == GeoClass::Tri {
                    2.0 / (ntriangle as f64)
                } else {
                    4.0 / (ntriangle as f64)
                };
                let mut pad = aux::gen_scratchpad_with_coords_aligned(kind);
                let nnode_total = kind.nnode() + kind.triangulate_extra_nnode();
                let mut xx = vec![0.0; nnode_total];
                let mut yy = vec![0.0; nnode_total];
                let zz = vec![0.0; nnode_total];
                let mut triangles = vec![vec![0; 3]; ntriangle];
                pad.triangulate(|t, i, m, x| {
                    xx[m] = x[0];
                    yy[m] = x[1];
                    triangles[t][i] = m;
                })
                .unwrap();
                for t in 0..ntriangle {
                    let a = triangles[t][0];
                    let b = triangles[t][1];
                    let c = triangles[t][2];
                    let xa = Point2d::from_slice(&[xx[a], yy[a]]);
                    let xb = Point2d::from_slice(&[xx[b], yy[b]]);
                    let xc = Point2d::from_slice(&[xx[c], yy[c]]);
                    let tri = Triangle2d::from_points(&xa, &xb, &xc);
                    approx_eq(tri.signed_area(), expected_area, 1e-14);
                    let (h0, h1, h2) = tri.internal_angles();
                    if kind.class() == GeoClass::Tri {
                        let mut angles = vec![h0, h1, h2];
                        angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
                        let (smallest, middle, largest) = (angles[0], angles[1], angles[2]);
                        approx_eq(smallest, forty_five, 1e-14);
                        approx_eq(middle, forty_five, 1e-14);
                        approx_eq(largest, ninety, 1e-14);
                    } else {
                        if t % 2 == 0 {
                            approx_eq(h0, ninety, 1e-14);
                            approx_eq(h1, forty_five, 1e-14);
                            approx_eq(h2, forty_five, 1e-14);
                        } else {
                            approx_eq(h0, forty_five, 1e-14);
                            approx_eq(h1, ninety, 1e-14);
                            approx_eq(h2, forty_five, 1e-14);
                        }
                    }
                }
                if SAVE_FIGURE {
                    let mut canvas = Canvas::new();
                    canvas.draw_triangles_3d(&xx, &yy, &zz, &triangles);
                    let mut plot = Plot::new();
                    pad.draw_shape(&mut plot, "", true, false).unwrap();
                    plot.add(&canvas);
                    plot.set_figure_size_points(600.0, 600.0)
                        .set_equal_axes(true)
                        .save(format!("/tmp/gemlab/test_pad_triangulation_works_{}.svg", kind.to_string()).as_str())
                        .unwrap();
                }
            }
        }
    }

    #[test]
    fn pad_triangulate_works_2() {
        let space_ndim = 3; // shells
        for kind in GeoKind::VALUES {
            let ntriangle = kind.triangulate_ntriangle();
            if ntriangle > 0 {
                let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);
                let nnode_total = kind.nnode() + kind.triangulate_extra_nnode();
                let mut xx = vec![0.0; nnode_total];
                let mut yy = vec![0.0; nnode_total];
                let mut zz = vec![0.0; nnode_total];
                let mut triangles = vec![vec![0; 3]; ntriangle];
                pad.triangulate(|t, i, m, x| {
                    xx[m] = x[0];
                    yy[m] = x[1];
                    zz[m] = x[2];
                    triangles[t][i] = m;
                })
                .unwrap();
                for t in 0..ntriangle {
                    let a = triangles[t][0];
                    let b = triangles[t][1];
                    let c = triangles[t][2];
                    let xa = Point2d::from_slice(&[xx[a], yy[a]]);
                    let xb = Point2d::from_slice(&[xx[b], yy[b]]);
                    let xc = Point2d::from_slice(&[xx[c], yy[c]]);
                    let tri = Triangle2d::from_points(&xa, &xb, &xc);
                    assert!(tri.signed_area() > 0.0);
                }
                if SAVE_FIGURE {
                    let mut canvas = Canvas::new();
                    canvas.draw_triangles_3d(&xx, &yy, &zz, &triangles);
                    let mut plot = Plot::new();
                    pad.draw_shape(&mut plot, "", true, false).unwrap();
                    plot.add(&canvas);
                    plot.set_figure_size_points(600.0, 600.0)
                        .set_equal_axes(true)
                        .save(format!("/tmp/gemlab/test_pad_triangulation_works_2_{}.svg", kind.to_string()).as_str())
                        .unwrap();
                }
            }
        }
    }
}

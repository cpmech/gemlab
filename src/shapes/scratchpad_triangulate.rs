#![allow(unused)]

use super::{GeoClass, Scratchpad};
use crate::StrError;
use russell_lab::{mat_vec_mul, vec_norm, Norm, Vector};

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
    use crate::shapes::scratchpad_testing::aux;
    use crate::shapes::GeoKind;
    use crate::shapes::Scratchpad;
    use plotpy::{Canvas, Plot};
    use russell_lab::math::{ONE_BY_3, SQRT_3};
    use russell_lab::{array_approx_eq, vec_approx_eq, Vector};

    const SAVE_FIGURE: bool = false;

    #[test]
    fn pad_triangulate_works() {
        let space_ndim = 3; // shells in 3D
                            // for kind in [GeoKind::Qua8] {
        for kind in GeoKind::VALUES {
            let ntriangle = kind.triangulate_ntriangle();
            if ntriangle > 0 {
                let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);
                // pad.xxt.set(2, 4, 3.0);
                // pad.xxt.set(2, 6, 5.0);
                println!("Testing kind {:?}", kind);
                let mut nnode_total = kind.nnode() + kind.triangulate_extra_nnode();
                let mut xx = vec![0.0; nnode_total];
                let mut yy = vec![0.0; nnode_total];
                let mut zz = vec![0.0; nnode_total];
                let mut triangles = vec![vec![0; 3]; ntriangle];
                pad.triangulate(|t, i, m, x| {
                    println!("t = {}, i = {}, m = {},  x = {:?}", t, i, m, x.as_data());
                    xx[m] = x[0];
                    yy[m] = x[1];
                    zz[m] = x[2];
                    triangles[t][i] = m;
                });
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
}

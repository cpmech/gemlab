#[cfg(test)]
pub mod aux {
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::util::PI;
    use crate::StrError;
    use plotpy::{Canvas, Curve, Plot};
    use russell_lab::generate2d;
    use russell_lab::Vector;

    pub const RMIN: f64 = 5.0;
    pub const RMAX: f64 = 10.0;
    pub const AMIN: f64 = 30.0 * PI / 180.0;
    pub const AMAX: f64 = 60.0 * PI / 180.0;
    pub const ZMIN: f64 = 0.0;
    pub const ZMAX: f64 = 4.0;

    /// Maps point coordinates
    ///
    /// * The shape is the area indicated with "?" or the edge with "%".
    /// * If class == Tri, the shape is the half (bottom) of the highlighted wedge.
    /// * In 3D, an extrusion is applied along the out-of-plane direction.
    ///
    /// ```text
    ///   |            /
    ///   |           / αmax
    ///   ***=---__  /
    ///   |         % _
    ///   |        % ? *._          ,
    ///   |       % ????? *.     ,-'
    ///   ***=-_ % ???????? *.,-' αmin
    ///   |     % - ?????? ,-'*
    ///   |    /    *.? ,-'    *
    ///   |   /      ,*'        *
    ///   |  /    ,-'  *         *
    ///   | /  ,-'      *         *
    ///   |/.-'         #         #
    ///   o ----------- # ------- # --> r
    ///               rmin       rmax
    /// ```
    ///
    /// Intermediary mapping:
    ///
    /// r(ξ₀,ξ₁,ξ₂) = rmin + (ξ₀ - ξ₀min) · Δr / Δξ₀
    /// α(ξ₀,ξ₁,ξ₂) = αmin + (ξ₁ - ξ₁min) · Δα / Δξ₁
    /// z(ξ₀,ξ₁,ξ₂) = ξ₂
    ///
    /// Cylindrical coordinates:
    ///
    /// x₀ := r · cos(α)
    /// x₁ := r · sin(α)
    /// x₂ := z
    pub fn map_point_coords(x: &mut Vector, ksi: &[f64], ksi_min: f64, ksi_del: f64) {
        assert_eq!(x.dim(), ksi.len());
        let r = RMIN + (ksi[0] - ksi_min) * (RMAX - RMIN) / ksi_del;
        let a = AMIN + (ksi[1] - ksi_min) * (AMAX - AMIN) / ksi_del;
        x[0] = r * f64::cos(a);
        x[1] = r * f64::sin(a);
        if x.dim() == 3 {
            x[2] = ZMIN + (ksi[2] - ksi_min) * (ZMAX - ZMIN) / ksi_del;
        }
    }

    /// Returns a new scratchpad with coordinates set
    ///
    /// Notes:
    ///
    /// * For cables, the line will be along AMAX
    /// * For shells, the surface will be at ZMAX
    pub fn gen_scratchpad_with_coords(space_ndim: usize, kind: GeoKind) -> Scratchpad {
        let geo_ndim = kind.ndim();
        let nnode = kind.nnode();
        let mut x = Vector::new(space_ndim);
        let mut ksi_aux = vec![0.0; space_ndim];
        let mut pad = Scratchpad::new(space_ndim, kind).unwrap();
        let (ksi_min, ksi_del) = kind.ksi_min_ksi_del();
        for m in 0..nnode {
            let ksi = kind.reference_coords(m);
            if geo_ndim == space_ndim {
                // SOLID case: 2D or 3D
                map_point_coords(&mut x, ksi, ksi_min, ksi_del);
            } else if geo_ndim == 1 && space_ndim == 2 {
                // CABLE case: line in 2D
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = 1.0;
                map_point_coords(&mut x, &ksi_aux, ksi_min, ksi_del);
            } else {
                // CABLE or SHELL in 3D
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = ksi[1];
                ksi_aux[2] = 1.0;
                map_point_coords(&mut x, &ksi_aux, ksi_min, ksi_del);
            }
            for j in 0..space_ndim {
                pad.set_xx(m, j, x[j])
            }
        }
        pad
    }

    /// Returns a new scratchpad with coordinates such that the shape has edges aligned to x-y-z
    ///
    /// **Important:** This function works with geo_ndim = 2 or 3 only (SOLID).
    ///
    /// Notes:
    ///
    /// * For triangles (tetrahedra), one edge (face) will cut the diagonal of the x-y(-z) space.
    /// * Qua and Hex will have real coordinates equal to the natural coordinates
    /// * Tri and Tet will be scaled using the natural coordinates
    /// * All edges parallel to the x,y,z axes will have lengths equal to 2.0
    /// * The cell will also be shifted such that the first vertex will be at (0.0,0.0,0.0)
    pub fn gen_scratchpad_with_coords_aligned(kind: GeoKind) -> Scratchpad {
        let geo_ndim = kind.ndim();
        assert!(geo_ndim > 1);
        let space_ndim = geo_ndim;
        let mut pad = Scratchpad::new(space_ndim, kind).unwrap();
        let (shift, scale) = if kind.is_tri_or_tet() { (0.0, 2.0) } else { (1.0, 1.0) };
        for m in 0..kind.nnode() {
            let ksi = kind.reference_coords(m);
            for j in 0..space_ndim {
                pad.set_xx(m, j, shift + scale * ksi[j])
            }
        }
        pad
    }

    /// Extracts edge 'e' with coordinates set
    ///
    /// # Panics
    ///
    /// A panic will occur if the shape does not have edges
    pub fn extract_edge(e: usize, pad: &Scratchpad) -> Result<Scratchpad, StrError> {
        let (space_ndim, geo_ndim) = pad.jacobian.dims();
        if geo_ndim == 1 {
            return Err("extract_edge requires geo_ndim to be 2 or 3");
        }
        let mut pad_edge = Scratchpad::new(space_ndim, pad.kind.edge_kind().unwrap())?;
        for i in 0..pad.kind.edge_nnode() {
            let m = pad.kind.edge_node_id(e, i);
            for j in 0..space_ndim {
                pad_edge.set_xx(i, j, pad.xxt[j][m]);
            }
        }
        Ok(pad_edge)
    }

    /// Extracts face 'f' with coordinates set
    ///
    /// # Panics
    ///
    /// A panic will occur if the shape does not have faces
    pub fn extract_face(f: usize, pad: &Scratchpad) -> Result<Scratchpad, StrError> {
        let (space_ndim, geo_ndim) = pad.jacobian.dims();
        if geo_ndim != 3 {
            return Err("extract_face requires geo_ndim to be 3");
        }
        let mut pad_face = Scratchpad::new(space_ndim, pad.kind.face_kind().unwrap())?;
        for i in 0..pad.kind.face_nnode() {
            let m = pad.kind.face_node_id(f, i);
            for j in 0..space_ndim {
                pad_face.set_xx(i, j, pad.xxt[j][m]);
            }
        }
        Ok(pad_face)
    }

    // ------ internal functions for testing ----------------------------------------------------------

    /// Generates the Canvas for the mapping used in tests
    #[cfg_attr(coverage_nightly, no_coverage)]
    pub fn _gen_canvas_mapping() -> Canvas {
        let mut canvas = Canvas::new();
        let color = "#bfbfbf";
        canvas
            .set_face_color("None")
            .set_edge_color(color)
            .set_line_width(2.0)
            .draw_circle(0.0, 0.0, RMIN);
        canvas.draw_circle(0.0, 0.0, RMAX);
        canvas.set_edge_color(color).set_line_width(2.0);
        canvas.draw_polyline(&[[0.0, 0.0], [RMAX * f64::cos(AMIN), RMAX * f64::sin(AMIN)]], false);
        canvas.draw_polyline(&[[0.0, 0.0], [RMAX * f64::cos(AMAX), RMAX * f64::sin(AMAX)]], false);
        canvas
    }

    /// Draws the points in the natural (right) and real (left) spaces
    #[cfg_attr(coverage_nightly, no_coverage)]
    pub fn _draw_point_coords_2d(ksi_min: f64, ksi_del: f64) -> Plot {
        const N: usize = 11;
        let mut natural_space = Curve::new();
        natural_space.set_line_style("None").set_marker_style("o");
        let mut real_space = Curve::new();
        real_space.set_line_style("None").set_marker_style("o");
        let ksi_max = ksi_min + ksi_del;
        let (ksi_0, ksi_1) = generate2d(ksi_min, ksi_max, ksi_min, ksi_max, N, N);
        let mut x = Vector::new(2);
        natural_space.points_begin();
        real_space.points_begin();
        for i in 0..N {
            for j in 0..N {
                natural_space.points_add(ksi_0[i][j], ksi_1[i][j]);
                map_point_coords(&mut x, &[ksi_0[i][j], ksi_1[i][j]], ksi_min, ksi_del);
                real_space.points_add(x[0], x[1]);
            }
        }
        natural_space.points_end();
        real_space.points_end();
        let mut plot = Plot::new();
        plot.set_subplot(1, 2, 1)
            .add(&real_space)
            .add(&_gen_canvas_mapping())
            .set_equal_axes(true)
            .set_range(-0.1, RMAX + 0.1, -0.1, RMAX + 0.1)
            .set_ticks_x(1.0, 0.0, "")
            .set_ticks_y(1.0, 0.0, "")
            .grid_and_labels("x", "y");
        plot.set_subplot(1, 2, 2)
            .add(&natural_space)
            .set_equal_axes(true)
            .set_labels("ξ0", "ξ1");
        plot
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::aux::{_draw_point_coords_2d, extract_edge, extract_face, gen_scratchpad_with_coords};
    use crate::shapes::op::testing::aux::{self, gen_scratchpad_with_coords_aligned};
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::StrError;
    use russell_chk::{assert_approx_eq, assert_deriv_approx_eq};

    // #[test]
    #[cfg_attr(coverage_nightly, no_coverage)]
    fn _draw_point_coords_2d_works() -> Result<(), StrError> {
        let mut plot = _draw_point_coords_2d(0.0, 1.0);
        plot.set_figure_size_points(1000.0, 600.0)
            .save("/tmp/gemlab/test_draw_point_coords_2d.svg")
    }

    #[test]
    fn gen_scratchpad_with_coords_works() {
        // SOLID in 2D
        let pad = gen_scratchpad_with_coords(2, GeoKind::Qua4);
        assert_approx_eq!(pad.xxt[0][0], aux::RMIN * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[1][0], aux::RMIN * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[0][1], aux::RMAX * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[1][1], aux::RMAX * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[0][2], aux::RMAX * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[1][2], aux::RMAX * f64::sin(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[0][3], aux::RMIN * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[1][3], aux::RMIN * f64::sin(aux::AMAX), 1e-15);

        // SOLID in 3D
        let pad = gen_scratchpad_with_coords(3, GeoKind::Hex8);
        assert_approx_eq!(pad.xxt[0][0], aux::RMIN * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[1][0], aux::RMIN * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[2][0], aux::ZMIN, 1e-15);
        assert_approx_eq!(pad.xxt[0][1], aux::RMAX * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[1][1], aux::RMAX * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[2][1], aux::ZMIN, 1e-15);
        assert_approx_eq!(pad.xxt[0][2], aux::RMAX * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[1][2], aux::RMAX * f64::sin(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[2][2], aux::ZMIN, 1e-15);
        assert_approx_eq!(pad.xxt[0][3], aux::RMIN * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[1][3], aux::RMIN * f64::sin(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[2][3], aux::ZMIN, 1e-15);

        assert_approx_eq!(pad.xxt[0][4], aux::RMIN * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[1][4], aux::RMIN * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[2][4], aux::ZMAX, 1e-15);
        assert_approx_eq!(pad.xxt[0][5], aux::RMAX * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[1][5], aux::RMAX * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[2][5], aux::ZMAX, 1e-15);
        assert_approx_eq!(pad.xxt[0][6], aux::RMAX * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[1][6], aux::RMAX * f64::sin(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[2][6], aux::ZMAX, 1e-15);
        assert_approx_eq!(pad.xxt[0][7], aux::RMIN * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[1][7], aux::RMIN * f64::sin(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[2][7], aux::ZMAX, 1e-15);

        // CABLE in 2D
        let pad = gen_scratchpad_with_coords(2, GeoKind::Lin2);
        assert_approx_eq!(pad.xxt[0][0], aux::RMIN * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[1][0], aux::RMIN * f64::sin(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[0][1], aux::RMAX * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[1][1], aux::RMAX * f64::sin(aux::AMAX), 1e-15);

        // SHELL in 3D
        let pad = gen_scratchpad_with_coords(3, GeoKind::Tri3);
        assert_approx_eq!(pad.xxt[0][0], aux::RMIN * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[1][0], aux::RMIN * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[2][0], aux::ZMAX, 1e-15);
        assert_approx_eq!(pad.xxt[0][1], aux::RMAX * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[1][1], aux::RMAX * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad.xxt[2][1], aux::ZMAX, 1e-15);
        assert_approx_eq!(pad.xxt[0][2], aux::RMIN * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[1][2], aux::RMIN * f64::sin(aux::AMAX), 1e-15);
        assert_approx_eq!(pad.xxt[2][2], aux::ZMAX, 1e-15);
    }

    #[test]
    fn gen_scratchpad_with_coords_aligned_works() {
        let pad = gen_scratchpad_with_coords_aligned(GeoKind::Tri3);
        assert_eq!(pad.xxt[0][0], 0.0);
        assert_eq!(pad.xxt[1][0], 0.0);
        assert_eq!(pad.xxt[0][1], 2.0);
        assert_eq!(pad.xxt[1][1], 0.0);
        assert_eq!(pad.xxt[0][2], 0.0);
        assert_eq!(pad.xxt[1][2], 2.0);

        let pad = gen_scratchpad_with_coords_aligned(GeoKind::Qua4);
        assert_eq!(pad.xxt[0][0], 0.0);
        assert_eq!(pad.xxt[1][0], 0.0);
        assert_eq!(pad.xxt[0][1], 2.0);
        assert_eq!(pad.xxt[1][1], 0.0);
        assert_eq!(pad.xxt[0][2], 2.0);
        assert_eq!(pad.xxt[1][2], 2.0);
        assert_eq!(pad.xxt[0][3], 0.0);
        assert_eq!(pad.xxt[1][3], 2.0);

        let pad = gen_scratchpad_with_coords_aligned(GeoKind::Hex8);
        assert_eq!(pad.xxt[0][0], 0.0);
        assert_eq!(pad.xxt[1][0], 0.0);
        assert_eq!(pad.xxt[2][0], 0.0);
        assert_eq!(pad.xxt[0][1], 2.0);
        assert_eq!(pad.xxt[1][1], 0.0);
        assert_eq!(pad.xxt[2][1], 0.0);
        assert_eq!(pad.xxt[0][2], 2.0);
        assert_eq!(pad.xxt[1][2], 2.0);
        assert_eq!(pad.xxt[2][2], 0.0);
        assert_eq!(pad.xxt[0][3], 0.0);
        assert_eq!(pad.xxt[1][3], 2.0);
        assert_eq!(pad.xxt[2][3], 0.0);

        assert_eq!(pad.xxt[0][4], 0.0);
        assert_eq!(pad.xxt[1][4], 0.0);
        assert_eq!(pad.xxt[2][4], 2.0);
        assert_eq!(pad.xxt[0][5], 2.0);
        assert_eq!(pad.xxt[1][5], 0.0);
        assert_eq!(pad.xxt[2][5], 2.0);
        assert_eq!(pad.xxt[0][6], 2.0);
        assert_eq!(pad.xxt[1][6], 2.0);
        assert_eq!(pad.xxt[2][6], 2.0);
        assert_eq!(pad.xxt[0][7], 0.0);
        assert_eq!(pad.xxt[1][7], 2.0);
        assert_eq!(pad.xxt[2][7], 2.0);
    }

    #[test]
    fn capture_some_errors() {
        let pad = Scratchpad::new(3, GeoKind::Lin2).unwrap();
        assert_eq!(
            extract_edge(0, &pad).err(),
            Some("extract_edge requires geo_ndim to be 2 or 3")
        );
        let pad = Scratchpad::new(3, GeoKind::Tri3).unwrap();
        assert_eq!(
            extract_face(0, &pad).err(),
            Some("extract_face requires geo_ndim to be 3")
        );
    }

    #[test]
    fn extract_edge_works() -> Result<(), StrError> {
        let pad = gen_scratchpad_with_coords(2, GeoKind::Qua4);
        let pad_edge = extract_edge(0, &pad)?;
        assert_eq!(pad_edge.kind, GeoKind::Lin2);
        assert_approx_eq!(pad_edge.xxt[0][0], aux::RMAX * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad_edge.xxt[1][0], aux::RMAX * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad_edge.xxt[0][1], aux::RMIN * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad_edge.xxt[1][1], aux::RMIN * f64::sin(aux::AMIN), 1e-15);
        Ok(())
    }

    #[test]
    fn extract_face_works() -> Result<(), StrError> {
        let pad = gen_scratchpad_with_coords(3, GeoKind::Hex8);
        let pad_face = extract_face(0, &pad)?;
        assert_eq!(pad_face.kind, GeoKind::Qua4);
        assert_approx_eq!(pad_face.xxt[0][0], aux::RMIN * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad_face.xxt[1][0], aux::RMIN * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad_face.xxt[2][0], aux::ZMIN, 1e-15);
        assert_approx_eq!(pad_face.xxt[0][1], aux::RMIN * f64::cos(aux::AMIN), 1e-15);
        assert_approx_eq!(pad_face.xxt[1][1], aux::RMIN * f64::sin(aux::AMIN), 1e-15);
        assert_approx_eq!(pad_face.xxt[2][1], aux::ZMAX, 1e-15);
        assert_approx_eq!(pad_face.xxt[0][2], aux::RMIN * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad_face.xxt[1][2], aux::RMIN * f64::sin(aux::AMAX), 1e-15);
        assert_approx_eq!(pad_face.xxt[2][2], aux::ZMAX, 1e-15);
        assert_approx_eq!(pad_face.xxt[0][3], aux::RMIN * f64::cos(aux::AMAX), 1e-15);
        assert_approx_eq!(pad_face.xxt[1][3], aux::RMIN * f64::sin(aux::AMAX), 1e-15);
        assert_approx_eq!(pad_face.xxt[2][3], aux::ZMIN, 1e-15);
        Ok(())
    }

    #[test]
    fn fn_interp_works() -> Result<(), StrError> {
        // loop over shapes
        for kind in GeoKind::VALUES {
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = gen_scratchpad_with_coords(space_ndim, kind);

            // loop over nodes
            const TOL: f64 = 1e-15;
            let nnode = kind.nnode();
            for m in 0..nnode {
                // get ξᵐ corresponding to node m
                let ksi = kind.reference_coords(m);

                // compute interpolation function Nⁿ(ξᵐ)
                (pad.fn_interp)(&mut pad.interp, ksi);

                // check: Nⁿ(ξᵐ) = 1 if m==n; 0 otherwise
                for n in 0..nnode {
                    if m == n {
                        assert_approx_eq!(pad.interp[n], 1.0, TOL);
                    } else {
                        assert_approx_eq!(pad.interp[n], 0.0, TOL);
                    }
                }
            }
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of N with respect to ξ => L (deriv) matrix
    struct ArgsNumDeriv {
        pad: Scratchpad,  // scratchpad
        at_ksi: Vec<f64>, // at reference coord value
        ksi: Vec<f64>,    // temporary reference coord
        m: usize,         // node index from 0 to nnode
        j: usize,         // dimension index from 0 to geom_ndim
    }

    // Computes Nᵐ(ξ) with variable v := ξⱼ
    fn nn_given_ksi(v: f64, args: &mut ArgsNumDeriv) -> f64 {
        args.ksi.copy_from_slice(&args.at_ksi);
        args.ksi[args.j] = v;
        (args.pad.fn_interp)(&mut args.pad.interp, &args.ksi);
        args.pad.interp[args.m]
    }

    #[test]
    fn calc_deriv_works() -> Result<(), StrError> {
        // kind and tolerances
        let problem = vec![
            // Lin
            (GeoKind::Lin2, 1e-13),
            (GeoKind::Lin3, 1e-13),
            (GeoKind::Lin4, 1e-10),
            (GeoKind::Lin5, 1e-10),
            // Tri
            (GeoKind::Tri3, 1e-12),
            (GeoKind::Tri6, 1e-12),
            (GeoKind::Tri10, 1e-10),
            (GeoKind::Tri15, 1e-9),
            // Qua
            (GeoKind::Qua4, 1e-13),
            (GeoKind::Qua8, 1e-12),
            (GeoKind::Qua9, 1e-13),
            (GeoKind::Qua12, 1e-10),
            (GeoKind::Qua16, 1e-10),
            (GeoKind::Qua17, 1e-10),
            // Tet
            (GeoKind::Tet4, 1e-12),
            (GeoKind::Tet10, 1e-12),
            (GeoKind::Tet20, 1e-10),
            // Hex
            (GeoKind::Hex8, 1e-13),
            (GeoKind::Hex20, 1e-12),
            (GeoKind::Hex32, 1e-10),
        ];

        // loop over shapes
        for (kind, tol) in problem {
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // set ξ within reference space
            let at_ksi = vec![0.25; geo_ndim];

            // compute all derivatives of interpolation functions with respect to ξ
            (pad.fn_deriv)(&mut pad.deriv, &at_ksi);

            // set arguments for numerical integration
            let args = &mut ArgsNumDeriv {
                pad: pad.clone(),
                at_ksi,
                ksi: vec![0.0; geo_ndim],
                m: 0,
                j: 0,
            };

            // check Lᵐ(ξ) = dNᵐ(ξ)/dξ
            for m in 0..kind.nnode() {
                args.m = m;
                for j in 0..geo_ndim {
                    args.j = j;
                    // Lᵐⱼ := dNᵐ/dξⱼ
                    assert_deriv_approx_eq!(pad.deriv[m][j], args.at_ksi[j], nn_given_ksi, args, tol);
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
pub mod aux {
    use crate::shapes::GeoKind;
    use crate::util::PI;
    use russell_lab::{Matrix, Vector};

    /// Maps point coordinates
    ///
    /// The shape is the area indicated with "?" or the edge with "%".
    /// If class == Tri, the shape is half of the highlighted wedge.
    /// In 3D, an extrusion is applied along the out-of-plane direction.
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
    /// r(ξ₀,ξ₁,ξ₂) = rmin + (ξ₀ - ξ₀min) * Δr / Δξ₀
    /// α(ξ₀,ξ₁,ξ₂) = αmin + (ξ₁ - ξ₁min) * Δα / Δξ₁
    /// z(ξ₀,ξ₁,ξ₂) = ξ₂
    ///
    /// Cylindrical coordinates:
    ///
    /// x₀ := r * cos(α)
    /// x₁ := r * sin(α)
    /// x₂ := z
    pub fn map_point_coords(x: &mut Vector, ksi: &[f64], kind: GeoKind) {
        const RMIN: f64 = 1.0;
        const RMAX: f64 = 10.0;
        const AMIN: f64 = 30.0 * PI / 180.0;
        const AMAX: f64 = 60.0 * PI / 180.0;
        assert_eq!(x.dim(), ksi.len());
        let (min_ksi, _, del_ksi) = kind.reference_limits();
        let r = RMIN + (ksi[0] - min_ksi) * (RMAX - RMIN) / del_ksi;
        let a = AMIN + (ksi[1] - min_ksi) * (AMAX - AMIN) / del_ksi;
        x[0] = r * f64::cos(a);
        x[1] = r * f64::sin(a);
        if x.dim() == 3 {
            x[2] = ksi[2];
        }
    }

    /// Generates a matrix of coordinates (transposed) according to the
    /// function `map_point_coords` above
    ///
    /// ```text
    ///      ┌                              ┐  superscript = node
    ///      | x⁰₀  x¹₀  x²₀  x³₀       xᴹ₀ |  subscript = dimension
    /// Xᵀ = | x⁰₁  x¹₁  x²₁  x³₁  ...  xᴹ₁ |
    ///      | x⁰₂  x¹₂  x²₂  x³₂       xᴹ₂ |
    ///      └                              ┘_(space_ndim,nnode)
    /// ```
    pub fn gen_coords_transp(space_ndim: usize, kind: GeoKind) -> Matrix {
        let geo_ndim = kind.ndim();
        let nnode = kind.nnode();
        let mut x = Vector::new(space_ndim);
        let mut ksi_aux = vec![0.0; space_ndim];
        let mut xxt = Matrix::new(space_ndim, nnode);
        for m in 0..nnode {
            let ksi = kind.reference_coords(m);
            if geo_ndim == space_ndim {
                map_point_coords(&mut x, ksi, kind);
            } else if geo_ndim == 1 && space_ndim == 2 {
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = 1.0;
                map_point_coords(&mut x, &ksi_aux, kind);
            } else {
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = ksi[1];
                ksi_aux[2] = 1.0;
                map_point_coords(&mut x, &ksi_aux, kind);
            }
            for j in 0..space_ndim {
                xxt[j][m] = x[j];
            }
        }
        xxt
    }
}

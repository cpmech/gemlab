use gemlab::integ::Gauss;
use gemlab::shapes::GeoClass;
use gemlab::StrError;
use russell_lab::{mat_pseudo_inverse, mat_to_mathematica, mat_to_static_array, Matrix};
use std::fmt::Write;
use std::fs::{self, File};
use std::io::Write as IoWrite;
use std::path::Path;

// Calculates the pseudo-inverse of the reference coordinates matrix
//
// This function calculates `PINV_HXI ← pinv(hat_ξ)`, where `HXI ← hat_ξ` is
// a matrix with the reference (natural) coordinates of the integration points,
// augmented by a column of ones. Here, `pinv` means the pseudo-inverse.
//
// # Output
//
// ```text
// TR_PINV_HXI ← transpose(PINV_HXI) ← transpose(pinv(hat_ξ)))
// ```
//
// This function returns a string representation of the **transpose** of `PINV_HXI`.
// The reason for returning the transpose is due to the way static arrays are
// handled in Rust. With the transpose, we can share a references to the static
// data from any GeoClass. To do so, we also save the transpose matrix with
// 4 columns, even if geo_ndim < 3. In this way, we can write:
//
// ```text
// let q: &'static [[f64; 4]] = match class {
//    ...
// }
// ```
//
// # Reference
//
// 1. Durand R and Farias MM (2014) A local extrapolation method for finite elements,
//    Advances in Engineering Software, 67:1-9 <https://doi.org/10.1016/j.advengsoft.2013.07.002>
fn calc_pseudo_inverse_ref_coords_mat(
    class: GeoClass,
    n_integ_point: usize,
    mathematica: bool, // for mathematica check
    digits: usize,     // for mathematica check
) -> String {
    // get type of integration points set
    let typ = match class {
        GeoClass::Lin => "LEGENDRE",
        GeoClass::Tri => {
            if n_integ_point == 6 || n_integ_point == 7 {
                "FELIPPA"
            } else {
                "INTERNAL"
            }
        }
        GeoClass::Qua => "LEGENDRE",
        GeoClass::Tet => {
            if n_integ_point > 5 {
                "FELIPPA"
            } else {
                "INTERNAL"
            }
        }
        GeoClass::Hex => {
            if n_integ_point == 6 || n_integ_point == 14 {
                "IRONS"
            } else {
                "LEGENDRE"
            }
        }
    };

    // get const name
    let name = match class {
        GeoClass::Lin => format!("TR_PINV_HXI_LIN_{}_{}", typ, n_integ_point),
        GeoClass::Tri => format!("TR_PINV_HXI_TRI_{}_{}", typ, n_integ_point),
        GeoClass::Qua => format!("TR_PINV_HXI_QUA_{}_{}", typ, n_integ_point),
        GeoClass::Tet => format!("TR_PINV_HXI_TET_{}_{}", typ, n_integ_point),
        GeoClass::Hex => format!("TR_PINV_HXI_HEX_{}_{}", typ, n_integ_point),
    };

    // integration points
    let gauss = Gauss::new_sized(class, n_integ_point).unwrap();

    // hat_ξ matrix (n_integ_point,geo_ndim+1) with the natural coordinates of the integration points, augmented by a column of ones.
    let geo_ndim = class.ndim();
    let mut hxi = Matrix::new(n_integ_point, geo_ndim + 1);
    for p in 0..n_integ_point {
        for d in 0..geo_ndim {
            hxi.set(p, d, gauss.coords(p)[d]);
        }
        hxi.set(p, geo_ndim, 1.0);
    }

    // print Mathematica code
    let mut buf = String::new();
    if mathematica {
        write!(&mut buf, "\n(* {} *)\n", name).unwrap();
        write!(&mut buf, "{}", mat_to_mathematica("hxi", &hxi)).unwrap();
    }

    // calculate pinv(hat_ξ) (geo_ndim+1,n_integ_point), the pseudo-inverse of hat_ξ
    let mut pinv_hxi = Matrix::new(geo_ndim + 1, n_integ_point);
    mat_pseudo_inverse(&mut pinv_hxi, &mut hxi).unwrap();

    // results
    if mathematica {
        // print Mathematica code
        write!(&mut buf, "{}", mat_to_mathematica("pinvHxi", &pinv_hxi)).unwrap();
        write!(
            &mut buf,
            "err = Max[Abs[pinvHxi - PseudoInverse[hxi]]];\n\
             Print[\"{}: err = \", err, \", check = \", err < 10^-{}]\n",
            name, digits,
        )
        .unwrap();
    } else {
        // generate transposed matrix with extra columns filled with zeros
        let mut tr_pinv_hxi = Matrix::new(n_integ_point, 4);
        for i in 0..n_integ_point {
            for j in 0..(geo_ndim + 1) {
                tr_pinv_hxi.set(i, j, pinv_hxi.get(j, i));
            }
        }

        // print Rust code
        write!(&mut buf, "\n{}", mat_to_static_array(&name, &tr_pinv_hxi)).unwrap();
    }
    buf
}

fn main() -> Result<(), StrError> {
    // for Mathematica, run with:
    //
    //    math < /tmp/gemlab/generate_gauss_extrap_data.m
    //
    let math = false; // mathematica
    let mut buf = String::new();

    // Lin
    if !math {
        buf.push_str("// -----------------------------------------------------------------------\n");
        buf.push_str("// -- LIN ----------------------------------------------------------------\n");
        buf.push_str("// -----------------------------------------------------------------------\n");
    }
    for n_integ_point in [1, 2, 3, 4, 5] {
        let res = calc_pseudo_inverse_ref_coords_mat(GeoClass::Lin, n_integ_point, math, 15);
        buf.push_str(&res);
    }

    // Tri
    if !math {
        buf.push_str("\n// -----------------------------------------------------------------------\n");
        buf.push_str("// -- TRI ----------------------------------------------------------------\n");
        buf.push_str("// -----------------------------------------------------------------------\n");
    }
    for n_integ_point in [1, 3, 4, 6, 7, 12, 16] {
        let res = calc_pseudo_inverse_ref_coords_mat(GeoClass::Tri, n_integ_point, math, 14);
        buf.push_str(&res);
    }

    // Qua
    if !math {
        buf.push_str("\n// -----------------------------------------------------------------------\n");
        buf.push_str("// -- QUA ----------------------------------------------------------------\n");
        buf.push_str("// -----------------------------------------------------------------------\n");
    }
    for n_integ_point in [1, 4, 9, 16] {
        let res = calc_pseudo_inverse_ref_coords_mat(GeoClass::Qua, n_integ_point, math, 15);
        buf.push_str(&res);
    }

    // Tet
    if !math {
        buf.push_str("\n// -----------------------------------------------------------------------\n");
        buf.push_str("// -- TET ----------------------------------------------------------------\n");
        buf.push_str("// -----------------------------------------------------------------------\n");
    }
    for n_integ_point in [1, 4, 5, 8, 14, 15, 24] {
        let res = calc_pseudo_inverse_ref_coords_mat(GeoClass::Tet, n_integ_point, math, 14);
        buf.push_str(&res);
    }

    // Hex
    if !math {
        buf.push_str("\n// -----------------------------------------------------------------------\n");
        buf.push_str("// -- HEX ----------------------------------------------------------------\n");
        buf.push_str("// -----------------------------------------------------------------------\n");
    }
    for n_integ_point in [6, 8, 14, 27, 64] {
        let res = calc_pseudo_inverse_ref_coords_mat(GeoClass::Hex, n_integ_point, math, 15);
        buf.push_str(&res);
    }

    // write file
    let ext = if math { "m" } else { "rs" };
    let name = format!("/tmp/gemlab/generate_gauss_extrap_data.{}", ext);
    let path = Path::new(&name);
    if let Some(p) = path.parent() {
        fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
    }
    let mut file = File::create(path).map_err(|_| "cannot create file")?;
    file.write_all(buf.as_bytes()).map_err(|_| "cannot write file")?;
    file.sync_all().map_err(|_| "cannot sync file")?;
    Ok(())
}

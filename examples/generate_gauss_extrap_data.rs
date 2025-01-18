use gemlab::integ::Gauss;
use gemlab::shapes::GeoClass;
use gemlab::StrError;
use russell_lab::{format_scientific, mat_pseudo_inverse, Matrix};
use std::fmt::Write;
use std::fs::{self, File};
use std::io::Write as IoWrite;
use std::path::Path;

fn write_mathematica_list(name: &str, a: &Matrix) -> String {
    let (nrow, ncol) = a.dims();
    let mut buf = String::new();
    write!(&mut buf, "{} = {{\n", name).unwrap();
    for i in 0..nrow {
        if i > 0 {
            write!(&mut buf, " }},\n").unwrap();
        }
        for j in 0..ncol {
            if j == 0 {
                write!(&mut buf, "{{ ").unwrap();
            } else {
                write!(&mut buf, ", ").unwrap();
            }
            let val = a.get(i, j);
            write!(&mut buf, "{:.17}", val).unwrap();
        }
    }
    write!(&mut buf, " }}\n").unwrap();
    write!(&mut buf, "}};\n").unwrap();
    buf
}

fn write_static_array(name: &str, a: &Matrix) -> String {
    let (nrow, ncol) = a.dims();
    let mut buf = String::new();
    write!(&mut buf, "#[rustfmt::skip]\n").unwrap();
    write!(&mut buf, "const {}: [[f64; {}]; {}] = [\n", name, ncol, nrow).unwrap();
    for i in 0..nrow {
        if i > 0 {
            write!(&mut buf, " ],\n").unwrap();
        }
        for j in 0..ncol {
            if j == 0 {
                write!(&mut buf, "\x20\x20\x20\x20[").unwrap();
            } else {
                write!(&mut buf, ",").unwrap();
            }
            let val = a.get(i, j);
            write!(&mut buf, "{}", format_scientific(val, 25, 17)).unwrap();
        }
    }
    write!(&mut buf, " ],\n").unwrap();
    write!(&mut buf, "];\n").unwrap();
    buf
}

// Calculates the pseudo-inverse of the reference (natural) coordinates matrix
//
// Returns the `ξ_hat_inv` matrix given by Eq. (32) of Reference #1.
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
    // set variable name
    let name = match class {
        GeoClass::Lin => format!("KSI_HAT_INV_LIN_{}", n_integ_point),
        GeoClass::Tri => format!("KSI_HAT_INV_TRI_{}", n_integ_point),
        GeoClass::Qua => format!("KSI_HAT_INV_QUA_{}", n_integ_point),
        GeoClass::Tet => format!("KSI_HAT_INV_TET_{}", n_integ_point),
        GeoClass::Hex => format!("KSI_HAT_INV_HEX_{}", n_integ_point),
    };

    // integration points
    let gauss = Gauss::new_sized(class, n_integ_point).unwrap();

    // ξ_hat matrix (n_integ_point,geo_ndim+1) with the natural coordinates of the integration points, augmented by a column of ones.
    // From Ref #1: ξ_hat is a matrix containing the local coordinates of the sampling (integration) points (Eq. (30) of Ref #1)
    let geo_ndim = class.ndim();
    let mut xh = Matrix::new(n_integ_point, geo_ndim + 1);
    for p in 0..n_integ_point {
        for d in 0..geo_ndim {
            xh.set(p, d, gauss.coords(p)[d]);
        }
        xh.set(p, geo_ndim, 1.0);
    }

    // print
    let mut buf = String::new();
    if mathematica {
        write!(&mut buf, "\n(* {} *)\n", name).unwrap();
        write!(&mut buf, "{}", write_mathematica_list("xh", &xh)).unwrap();
    }

    // calculate ξ_hat_inv (geo_ndim+1,n_integ_point), the pseudo-inverse of ξ_hat
    let mut xhi = Matrix::new(geo_ndim + 1, n_integ_point);
    mat_pseudo_inverse(&mut xhi, &mut xh).unwrap();

    // print
    if mathematica {
        write!(&mut buf, "{}", write_mathematica_list("xhi", &xhi)).unwrap();
        write!(
            &mut buf,
            "err = Max[Abs[xhi - PseudoInverse[xh]]];\n\
             Print[\"{}: err = \", err, \", check = \", err < 10^-{}]\n",
            name, digits,
        )
        .unwrap();
    } else {
        write!(&mut buf, "{}", write_static_array(&name, &xhi)).unwrap();
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
        buf.push_str("// -----------------------------------------------------------------------\n\n");
    }
    for n_integ_point in [1, 2, 3, 4, 5] {
        let res = calc_pseudo_inverse_ref_coords_mat(GeoClass::Lin, n_integ_point, math, 15);
        buf.push_str(&res);
    }

    // Tri
    if !math {
        buf.push_str("\n// -----------------------------------------------------------------------\n");
        buf.push_str("// -- TRI ----------------------------------------------------------------\n");
        buf.push_str("// -----------------------------------------------------------------------\n\n");
    }
    for n_integ_point in [1, 3, 4, 6, 7, 12, 16] {
        let res = calc_pseudo_inverse_ref_coords_mat(GeoClass::Tri, n_integ_point, math, 14);
        buf.push_str(&res);
    }

    // Qua
    if !math {
        buf.push_str("\n// -----------------------------------------------------------------------\n");
        buf.push_str("// -- QUA ----------------------------------------------------------------\n");
        buf.push_str("// -----------------------------------------------------------------------\n\n");
    }
    for n_integ_point in [1, 4, 9, 16] {
        let res = calc_pseudo_inverse_ref_coords_mat(GeoClass::Qua, n_integ_point, math, 15);
        buf.push_str(&res);
    }

    // Tet
    if !math {
        buf.push_str("\n// -----------------------------------------------------------------------\n");
        buf.push_str("// -- TET ----------------------------------------------------------------\n");
        buf.push_str("// -----------------------------------------------------------------------\n\n");
    }
    for n_integ_point in [1, 4, 5, 8, 14, 15, 24] {
        let res = calc_pseudo_inverse_ref_coords_mat(GeoClass::Tet, n_integ_point, math, 14);
        buf.push_str(&res);
    }

    // Hex
    if !math {
        buf.push_str("\n// -----------------------------------------------------------------------\n");
        buf.push_str("// -- HEX ----------------------------------------------------------------\n");
        buf.push_str("// -----------------------------------------------------------------------\n\n");
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

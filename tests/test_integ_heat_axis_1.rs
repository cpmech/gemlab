use gemlab::integ::{self, AnalyticalTri3, CommonArgs};
use gemlab::prelude::*;
use gemlab::shapes::{GeoClass, GeoKind, Scratchpad, Tri3};
use russell_lab::{mat_approx_eq, vec_approx_eq, Matrix, Vector};

// Example 5.6.1 from Lewis et al. Page 146
// Lewis RW, Nithiarasu R, Seetharamu KN (2004) Fundamentals of
// the Finite Element Method for Heat and Fluid Flow, Wiley.

// Returns pad_tri, pad_conv, and pad_flux
fn generate_pads() -> (Scratchpad, Scratchpad, Scratchpad, AnalyticalTri3) {
    #[rustfmt::skip]
    let mesh = Mesh {
        ndim: 2,
        points: vec![
            Point { id: 0, marker: 0, coords: vec![15.0, 10.0] },
            Point { id: 1, marker: 0, coords: vec![25.0, 10.0] },
            Point { id: 2, marker: 0, coords: vec![20.0, 12.0] },
        ],
        cells: vec![
            Cell { id: 0, attribute: 1, kind: GeoKind::Tri3, points: vec![0, 1, 2] },
        ],
    };
    let mut pad_tri = Scratchpad::new(2, GeoKind::Tri3).unwrap();
    let mut pad_conv = Scratchpad::new(2, GeoKind::Lin2).unwrap();
    let mut pad_flux = Scratchpad::new(2, GeoKind::Lin2).unwrap();
    mesh.set_pad(&mut pad_tri, &mesh.cells[0].points);
    mesh.set_pad(&mut pad_conv, &Tri3::EDGE_NODE_IDS[0]);
    mesh.set_pad(&mut pad_flux, &Tri3::EDGE_NODE_IDS[1]);
    let ana = AnalyticalTri3::new(&pad_tri);
    (pad_tri, pad_conv, pad_flux, ana)
}

#[test]
fn test_integ_heat_axis_1() {
    // generate scratchpads
    let (mut pad_tri, mut pad_conv, mut pad_flux, ana) = generate_pads();
    println!("area = {}", ana.area);
    println!("tri: coords =\n{}", pad_tri.xxt);
    println!("conv: coords =\n{}", pad_conv.xxt);
    println!("flux: coords =\n{}", pad_flux.xxt);

    // tri: arguments and integration points
    let ips_tri = integ::points(GeoClass::Tri, 3).unwrap(); // One integration point is enough
    let mut args_tri = CommonArgs::new(&mut pad_tri, ips_tri);
    args_tri.axisymmetric = true;

    // lin: sides of triangle
    let ips_lin = integ::points(GeoClass::Lin, 2).unwrap(); // IMPORTANT: we need at least 2 integration points
    let mut args_conv = CommonArgs::new(&mut pad_conv, ips_lin);
    let mut args_flux = CommonArgs::new(&mut pad_flux, ips_lin);
    args_conv.axisymmetric = true;
    args_flux.axisymmetric = true;

    // conductivity matrix
    let mut kk = Matrix::new(3, 3);
    let (kx, ky) = (2.0, 2.0);
    let kk_correct = ana.mat_03_btb(kx, ky, true);
    integ::mat_03_btb(&mut kk, &mut args_tri, |tt, _, _, _| {
        tt.sym_set(0, 0, kx);
        tt.sym_set(1, 1, ky);
        Ok(())
    })
    .unwrap();
    println!("kk = \n{:.1}", kk);
    println!("kk_correct = \n{:.1}", kk_correct);
    vec_approx_eq(kk.as_data(), kk_correct.as_data(), 1e-13);

    // convection on side 0
    let h_conv = 1.2;
    let mut kc = Matrix::new(2, 2);
    let kc_correct = ana.mat_01_nsn_bry(0, h_conv, true);
    integ::mat_01_nsn_bry(&mut kc, &mut args_conv, |_, _, _| Ok(h_conv)).unwrap();
    println!("kc = \n{}", kc);
    println!("kc_correct = \n{}", kc_correct);
    mat_approx_eq(&kc, &kc_correct, 1e-13);

    // sum conductivity and convection components
    for i in 0..2 {
        let m = GeoKind::Tri3.edge_node_id(0, i); // local(edge)-to-local(element)
        for j in 0..2 {
            let n = GeoKind::Tri3.edge_node_id(0, j); // local(edge)-to-local(element)
            kk.add(m, n, kc.get(i, j));
        }
    }
    let kk_lewis = Matrix::from(&[[99.0, 61.0, -50.0], [61.0, 119.0, -50.0], [-50.0, -50.0, 100.0]]);
    println!("kk_final = \n{:.1}", kk);
    println!("kk_lewis = \n{}", kk_lewis);
    mat_approx_eq(&kk, &kk_lewis, 1e-13);

    // source vector
    let source = 1.2;
    let mut f = Vector::new(3);
    integ::vec_01_ns(&mut f, &mut args_tri, |_, _| Ok(source)).unwrap();
    let f_source_correct = ana.vec_01_ns(source, true);
    println!("f_source = \n{}", f);
    println!("f_source_correct = \n{}", f_source_correct);
    vec_approx_eq(f.as_data(), f_source_correct.as_data(), 1e-13);

    // convection on side 0
    let mut f_conv = Vector::new(2);
    let temp_env = 30.0;
    integ::vec_01_ns(&mut f_conv, &mut args_conv, |_, _| Ok(h_conv * temp_env)).unwrap();
    let f_conv_correct = ana.vec_01_ns_bry(0, h_conv * temp_env, true);
    println!("f_conv = \n{}", f_conv);
    println!("f_conv_correct = \n{}", f_conv_correct);
    vec_approx_eq(f_conv.as_data(), f_conv_correct.as_data(), 1e-13);

    // flux on side 1
    let mut f_flux = Vector::new(2);
    let flux = 1.0;
    integ::vec_01_ns(&mut f_flux, &mut args_flux, |_, _| Ok(flux)).unwrap();
    let f_flux_correct = ana.vec_01_ns_bry(1, flux, true);
    println!("f_flux = \n{}", f_flux);
    println!("f_flux_correct = \n{}", f_flux_correct);
    vec_approx_eq(f_flux.as_data(), f_flux_correct.as_data(), 1e-13);

    // f-total
    for i in 0..2 {
        let m = GeoKind::Tri3.edge_node_id(0, i); // local(edge)-to-local(element)
        f[m] += f_conv[i];
        let m = GeoKind::Tri3.edge_node_id(1, i); // local(edge)-to-local(element)
        f[m] += f_flux[i];
    }
    println!("f_total = \n{}", f);

    let (ri, rj, rk) = (15.0, 25.0, 20.0);
    let (yi, yj, yk) = (10.0, 10.0, 12.0);
    let (lij, ljk) = (
        f64::sqrt((ri - rj) * (ri - rj) + (yi - yj) * (yi - yj)),
        f64::sqrt((rj - rk) * (rj - rk) + (yj - yk) * (yj - yk)),
    );
    println!("lij = {:?}, ljk = {:?}", lij, ljk);
    let (aa, gg, h, ta, q) = (ana.area, source, h_conv, temp_env, flux);
    println!("A = {}, G = {}, h = {}, Ta = {}, q = {}", aa, gg, h, ta, q);
    #[rustfmt::skip]
    let f_lewis= Vector::from(&[
        (gg*aa/12.0)*(2.0*ri+rj+rk) + 0.0                     + (h*ta*lij/6.0)*(2.0*ri+rj),
        (gg*aa/12.0)*(ri+2.0*rj+rk) + (q*ljk/6.0)*(2.0*rj+rk) + (h*ta*lij/6.0)*(ri+2.0*rj),
        (gg*aa/12.0)*(ri+rj+2.0*rk) + (q*ljk/6.0)*(rj+2.0*rk) + 0.0,
    ]);
    println!("f_lewis = \n{}", f_lewis);
    vec_approx_eq(f.as_data(), f_lewis.as_data(), 1e-13);
}

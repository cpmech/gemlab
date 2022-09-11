use gemlab::integ::{self, AnalyticalTri3, CommonArgs};
use gemlab::shapes::{GeoKind, Scratchpad};
use russell_chk::vec_approx_eq;
use russell_lab::{mat_approx_eq, Matrix};

/// Returns the scratchpad for a Tri3's side
///
/// **Important:** `side` must be 0, 1, or 2
fn tri3_extract_pad_side(pad: &Scratchpad, side: usize) -> Scratchpad {
    let mut pad_side = Scratchpad::new(2, GeoKind::Lin2).unwrap();
    for i in 0..2 {
        let m = pad.kind.edge_node_id(side, i);
        pad_side.set_xx(i, 0, pad.xxt[0][m]);
        pad_side.set_xx(i, 1, pad.xxt[1][m]);
    }
    pad_side
}

#[test]
fn test_integ_heat_eqs() {
    // Example 5.6.1 from Lewis et al. Page 146
    // Lewis RW, Nithiarasu R, Seetharamu KN (2004) Fundamentals of
    // the Finite Element Method for Heat and Fluid Flow, Wiley.
    //
    // Note: Convection at the bottom side of the triangle
    let mut pad = Scratchpad::new(2, GeoKind::Tri3).unwrap();
    pad.set_xx(0, 0, 15.0);
    pad.set_xx(0, 1, 10.0);
    pad.set_xx(1, 0, 25.0);
    pad.set_xx(1, 1, 10.0);
    pad.set_xx(2, 0, 20.0);
    pad.set_xx(2, 1, 12.0);

    // conductivity matrix
    let mut kk = Matrix::new(3, 3);
    let ana = AnalyticalTri3::new(&pad);
    let (kx, ky) = (2.0, 2.0);
    let kk_correct = ana.mat_03_btb(kx, ky, true);
    let class = pad.kind.class();
    let ips = integ::points(class, 1).unwrap();
    let mut args = CommonArgs::new(&mut pad, ips);
    args.axisymmetric = true;
    integ::mat_03_btb(&mut kk, &mut args, |tt, _, _, _| {
        tt.sym_set(0, 0, kx);
        tt.sym_set(1, 1, ky);
        Ok(())
    })
    .unwrap();
    println!("kk = \n{:.1}", kk);
    println!("kk_correct = \n{:.1}", kk_correct);
    vec_approx_eq(kk.as_data(), kk_correct.as_data(), 1e-14);

    // convection matrix
    let side = 0;
    let h_convection = 1.2;
    let mut kc = Matrix::new(2, 2);
    let mut pad_side = tri3_extract_pad_side(&pad, side);
    let kc_correct = ana.mat_01_nsn_bry(0, h_convection, true);
    let class = pad_side.kind.class();
    let ips = integ::points(class, 2).unwrap(); // IMPORTANT: we need at least 2 integration points here
    let mut args_side = CommonArgs::new(&mut pad_side, ips);
    args_side.axisymmetric = true;
    integ::mat_01_nsn_bry(&mut kc, &mut args_side, |_, _, _| Ok(h_convection)).unwrap();
    println!("kc = \n{}", kc);
    println!("kc_correct = \n{}", kc_correct);
    mat_approx_eq(&kc, &kc_correct, 1e-15);

    // add conductivity and convection components
    for i in 0..2 {
        let m = pad.kind.edge_node_id(side, i); // local(edge)-to-local(element)
        for j in 0..2 {
            let n = pad.kind.edge_node_id(side, j); // local(edge)-to-local(element)
            kk[m][n] += kc[i][j];
        }
    }
    let kk_lewis = Matrix::from(&[[99.0, 61.0, -50.0], [61.0, 119.0, -50.0], [-50.0, -50.0, 100.0]]);
    println!("kk_final = \n{}", kk);
    println!("kk_lewis = \n{}", kk_lewis);
    mat_approx_eq(&kk, &kk_lewis, 1e-15);
}

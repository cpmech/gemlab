use gemlab::integ;
use gemlab::shapes::{GeoKind, Scratchpad};
use gemlab::StrError;
use russell_chk::vec_approx_eq;
use russell_lab::Vector;
use russell_tensor::{Mandel, Tensor2};

fn main() -> Result<(), StrError> {
    // shape and state
    let space_ndim = 2;
    let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri3).unwrap();
    pad.set_xx(0, 0, 0.0);
    pad.set_xx(0, 1, 0.0);
    pad.set_xx(1, 0, 3.0);
    pad.set_xx(1, 1, 0.0);
    pad.set_xx(2, 0, 0.0);
    pad.set_xx(2, 1, 4.0);

    // analytical solutions
    let ana = integ::AnalyticalTri3::new(&pad);

    // shape times scalar, returns vector 'a'
    //
    //      ⌠    → →     →
    // aᵐ = │ Nᵐ(x(ξ)) s(x) dΩ
    //      ⌡
    //      Ωₑ
    //
    // with s(x) = 18
    //
    // solution (A = 3·4/2 = 6):
    //
    //          ┌   ┐   ┌    ┐
    //     18 A │ 1 │   │ 36 │
    // a = ———— │ 1 │ = │ 36 │
    //       3  │ 1 │   │ 36 │
    //          └   ┘   └    ┘
    let nnode = pad.kind.nnode();
    let ips = integ::default_points(pad.kind);
    let mut a = Vector::filled(nnode, 0.0);
    let mut args = integ::CommonArgs::new(&mut pad, ips);
    integ::vec_01_ns(&mut a, &mut args, |_, _| Ok(18.0))?;
    assert_eq!(
        format!("{:.1}", a),
        "┌      ┐\n\
         │ 36.0 │\n\
         │ 36.0 │\n\
         │ 36.0 │\n\
         └      ┘"
    );

    // check
    let a_correct = ana.vec_01_ns(18.0, false);
    vec_approx_eq(a.as_data(), a_correct.as_data(), 1e-14);

    // shape times vector, returns vector 'b'
    //
    // →    ⌠    → →   → →
    // bᵐ = │ Nᵐ(x(ξ)) v(x) dΩ
    //      ⌡
    //      Ωₑ
    //
    // with v(x) = {12, 12}
    //
    // solution (A = 3·4/2 = 6):
    //          ┌   ┐   ┌    ┐
    //          │ 1 │   │ 24 │
    //          │ 1 │   │ 24 │
    //     12 A │ 1 │   │ 24 │
    // b = ———— │ 1 │ = │ 24 │
    //       3  │ 1 │   │ 24 │
    //          │ 1 │   │ 24 │
    //          └   ┘   └    ┘
    // ```
    let mut b = Vector::filled(nnode * space_ndim, 0.0);
    integ::vec_02_nv(&mut b, &mut args, |v, _, _| {
        v[0] = 12.0;
        v[1] = 12.0;
        Ok(())
    })?;
    assert_eq!(
        format!("{}", b),
        "┌    ┐\n\
         │ 24 │\n\
         │ 24 │\n\
         │ 24 │\n\
         │ 24 │\n\
         │ 24 │\n\
         │ 24 │\n\
         └    ┘"
    );

    // check
    let b_correct = ana.vec_02_nv(12.0, 12.0);
    vec_approx_eq(b.as_data(), b_correct.as_data(), 1e-14);

    // vector dot gradient, returns vector 'c'
    //
    //      ⌠ → →    →  → →
    // cᵐ = │ w(x) · Bᵐ(x(ξ)) dΩ
    //      ⌡
    //      Ωₑ
    //
    // with w(x) = {-2.0, 4.0}
    //
    // solution:
    //     ┌    ┐
    //     │ -2 │
    // c = │ -4 │
    //     │  6 │
    //     └    ┘
    let mut c = Vector::filled(nnode, 0.0);
    integ::vec_03_vb(&mut c, &mut args, |w, _, _, _| {
        w[0] = -2.0;
        w[1] = 4.0;
        Ok(())
    })?;
    assert_eq!(
        format!("{}", c),
        "┌    ┐\n\
         │ -2 │\n\
         │ -4 │\n\
         │  6 │\n\
         └    ┘"
    );

    // check
    let c_correct = ana.vec_03_vb(-2.0, 4.0);
    vec_approx_eq(c.as_data(), c_correct.as_data(), 1e-15);

    // tensor dot gradient, returns vector 'd'
    //
    // →    ⌠   →    →  → →
    // dᵐ = │ σ(x) · Bᵐ(x(ξ)) dΩ
    //      ⌡ ▔
    //      Ωₑ
    //
    // with σ(x) = {σ₀₀, σ₁₁, σ₀₁√2}
    // and σ₀₀=3.0, σ₁₁=2.0, σ₀₁=0.5
    // σ₂₂ is ignored
    //
    // solution:
    //     ┌     ┐
    //     │ -15 │
    //     │ -10 │
    //     │  12 │
    // d = │   4 │
    //     │   3 │
    //     │   6 │
    //     └     ┘
    let (s00, s11, s01) = (6.0, 4.0, 2.0);
    let mut d = Vector::filled(nnode * space_ndim, 0.0);
    integ::vec_04_tb(&mut d, &mut args, |sig, _, _, _| {
        sig.sym_set(0, 0, s00);
        sig.sym_set(1, 1, s11);
        sig.sym_set(0, 1, s01);
        Ok(())
    })?;
    assert_eq!(
        format!("{}", d),
        "┌     ┐\n\
         │ -15 │\n\
         │ -10 │\n\
         │  12 │\n\
         │   4 │\n\
         │   3 │\n\
         │   6 │\n\
         └     ┘"
    );

    // check
    let sig = Tensor2::from_matrix(
        &[[s00, s01, 0.0], [s01, s11, 0.0], [0.0, 0.0, 0.0]],
        Mandel::Symmetric2D,
    )
    .unwrap();
    let d_correct = ana.vec_04_tb(&sig, false);
    vec_approx_eq(d.as_data(), d_correct.as_data(), 1e-15);
    Ok(())
}

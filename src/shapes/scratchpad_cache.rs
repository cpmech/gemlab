use russell_lab::{Matrix, Vector};

/// Maximum number of Gauss/reference points cached per [`super::Scratchpad`]
///
/// Matches the largest integration rule in `crate::integ::Gauss` (`Hex32`, 64 points) so every
/// standard element loop hits on every pass after the first while bounding memory on arbitrary
/// `ksi` sweeps (such as `Scratchpad::approximate_ksi`).
pub(crate) const MAX_CACHE_POINTS: usize = 64;

/// Converts up to `geo_ndim` reference coordinates `ksi` into a bitwise lookup key
#[inline]
pub(crate) fn ksi_key(geo_ndim: usize, ksi: &[f64]) -> [u64; 3] {
    match geo_ndim {
        1 => [ksi[0].to_bits(), 0, 0],
        2 => [ksi[0].to_bits(), ksi[1].to_bits(), 0],
        _ => [ksi[0].to_bits(), ksi[1].to_bits(), ksi[2].to_bits()],
    }
}

/// Bounded per-pad cache of `calc_jacobian` and `calc_gradient` evaluations keyed by `ksi`
///
/// Entries are laid out in flat contiguous buffers indexed by slot (`0..len <= MAX_CACHE_POINTS`)
/// so a cache hit requires one key comparison at a sequential cursor followed by `copy_from_slice`
/// into the `Scratchpad` arrays, with zero per-point heap allocations.
pub(crate) struct PadCache {
    len: usize,
    cursor: usize,
    nnode: usize,
    geo_ndim: usize,
    space_ndim: usize,
    geo_nnode: usize,
    mixed: bool,
    keys: Vec<[u64; 3]>,
    has_grad: Vec<bool>,
    det_jac: Vec<f64>,
    interp: Vec<f64>,
    deriv: Vec<f64>,
    deriv_geo: Vec<f64>,
    jacobian: Vec<f64>,
    inv_jacobian: Vec<f64>,
    gradient: Vec<f64>,
}

impl PadCache {
    /// Creates an empty cache sized for the geometry and field dimensions of a `Scratchpad`
    pub(crate) fn new(space_ndim: usize, geo_ndim: usize, nnode: usize, geo_nnode: usize, mixed: bool) -> Self {
        Self {
            len: 0,
            cursor: 0,
            nnode,
            geo_ndim,
            space_ndim,
            geo_nnode,
            mixed,
            keys: Vec::new(),
            has_grad: Vec::new(),
            det_jac: Vec::new(),
            interp: Vec::new(),
            deriv: Vec::new(),
            deriv_geo: Vec::new(),
            jacobian: Vec::new(),
            inv_jacobian: Vec::new(),
            gradient: Vec::new(),
        }
    }

    /// Discards all cached entries while retaining buffer capacities
    #[inline]
    pub(crate) fn clear(&mut self) {
        self.len = 0;
        self.cursor = 0;
    }

    /// Looks up a cached slot for `key`, checking the sequential cursor first
    #[inline]
    pub(crate) fn find(&mut self, key: [u64; 3]) -> Option<usize> {
        if self.len == 0 {
            return None;
        }
        if self.keys[self.cursor] == key {
            let idx = self.cursor;
            self.cursor = (idx + 1) % self.len;
            return Some(idx);
        }
        for idx in 0..self.len {
            if self.keys[idx] == key {
                self.cursor = (idx + 1) % self.len;
                return Some(idx);
            }
        }
        None
    }

    /// Returns whether `slot` already holds full `calc_gradient` outputs
    #[inline]
    pub(crate) fn has_gradient(&self, slot: usize) -> bool {
        self.has_grad[slot]
    }

    /// Allocates or evicts a slot for `key` and marks its gradient data absent
    pub(crate) fn alloc_slot(&mut self, key: [u64; 3]) -> usize {
        if self.len < MAX_CACHE_POINTS {
            let idx = self.len;
            self.len += 1;
            if self.len > self.keys.len() {
                let deriv_stride = self.nnode * self.geo_ndim;
                let deriv_geo_stride = if self.mixed { self.geo_nnode * self.geo_ndim } else { 0 };
                let jac_stride = self.space_ndim * self.geo_ndim;
                let inv_stride = if self.geo_ndim == self.space_ndim {
                    self.space_ndim * self.space_ndim
                } else {
                    0
                };
                let grad_stride = if self.geo_ndim == self.space_ndim {
                    self.nnode * self.space_ndim
                } else {
                    0
                };
                self.keys.resize(self.len, [0; 3]);
                self.has_grad.resize(self.len, false);
                self.det_jac.resize(self.len, 0.0);
                self.interp.resize(self.len * self.nnode, 0.0);
                self.deriv.resize(self.len * deriv_stride, 0.0);
                self.deriv_geo.resize(self.len * deriv_geo_stride, 0.0);
                self.jacobian.resize(self.len * jac_stride, 0.0);
                self.inv_jacobian.resize(self.len * inv_stride, 0.0);
                self.gradient.resize(self.len * grad_stride, 0.0);
            }
            self.keys[idx] = key;
            self.has_grad[idx] = false;
            self.cursor = 0;
            idx
        } else {
            let idx = self.cursor % MAX_CACHE_POINTS;
            self.keys[idx] = key;
            self.has_grad[idx] = false;
            self.cursor = (idx + 1) % MAX_CACHE_POINTS;
            idx
        }
    }

    /// Stores `calc_jacobian` outputs in `slot`
    pub(crate) fn store_jacobian(
        &mut self,
        slot: usize,
        det_jac: f64,
        deriv: &Matrix,
        deriv_geo: &Matrix,
        jacobian: &Matrix,
        inv_jacobian: &Matrix,
    ) {
        self.det_jac[slot] = det_jac;
        let jac_stride = self.space_ndim * self.geo_ndim;
        self.jacobian[slot * jac_stride..(slot + 1) * jac_stride].copy_from_slice(jacobian.as_data());

        if self.mixed {
            let stride = self.geo_nnode * self.geo_ndim;
            self.deriv_geo[slot * stride..(slot + 1) * stride].copy_from_slice(deriv_geo.as_data());
        } else {
            let stride = self.nnode * self.geo_ndim;
            self.deriv[slot * stride..(slot + 1) * stride].copy_from_slice(deriv.as_data());
        }

        if self.geo_ndim == self.space_ndim {
            let inv_stride = self.space_ndim * self.space_ndim;
            self.inv_jacobian[slot * inv_stride..(slot + 1) * inv_stride].copy_from_slice(inv_jacobian.as_data());
        }
    }

    /// Restores `calc_jacobian` outputs from `slot` without modifying `interp` or `gradient`
    #[inline]
    pub(crate) fn load_jacobian(
        &self,
        slot: usize,
        deriv: &mut Matrix,
        deriv_geo: &mut Matrix,
        jacobian: &mut Matrix,
        inv_jacobian: &mut Matrix,
    ) -> f64 {
        let jac_stride = self.space_ndim * self.geo_ndim;
        jacobian
            .as_mut_data()
            .copy_from_slice(&self.jacobian[slot * jac_stride..(slot + 1) * jac_stride]);

        if self.mixed {
            let stride = self.geo_nnode * self.geo_ndim;
            deriv_geo
                .as_mut_data()
                .copy_from_slice(&self.deriv_geo[slot * stride..(slot + 1) * stride]);
        } else {
            let stride = self.nnode * self.geo_ndim;
            deriv
                .as_mut_data()
                .copy_from_slice(&self.deriv[slot * stride..(slot + 1) * stride]);
        }

        if self.geo_ndim == self.space_ndim {
            let inv_stride = self.space_ndim * self.space_ndim;
            inv_jacobian
                .as_mut_data()
                .copy_from_slice(&self.inv_jacobian[slot * inv_stride..(slot + 1) * inv_stride]);
        }

        self.det_jac[slot]
    }

    /// Stores the remaining `calc_gradient` outputs in `slot` and marks `has_grad[slot] = true`
    pub(crate) fn store_gradient_part(&mut self, slot: usize, interp: &Vector, deriv: &Matrix, gradient: &Matrix) {
        self.interp[slot * self.nnode..(slot + 1) * self.nnode].copy_from_slice(interp.as_data());
        if self.mixed {
            let stride = self.nnode * self.geo_ndim;
            self.deriv[slot * stride..(slot + 1) * stride].copy_from_slice(deriv.as_data());
        }
        let grad_stride = self.nnode * self.space_ndim;
        self.gradient[slot * grad_stride..(slot + 1) * grad_stride].copy_from_slice(gradient.as_data());
        self.has_grad[slot] = true;
    }

    /// Restores full `calc_gradient` outputs from `slot` and returns `det_jac`
    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn load_gradient(
        &self,
        slot: usize,
        interp: &mut Vector,
        deriv: &mut Matrix,
        deriv_geo: &mut Matrix,
        jacobian: &mut Matrix,
        inv_jacobian: &mut Matrix,
        gradient: &mut Matrix,
    ) -> f64 {
        interp
            .as_mut_data()
            .copy_from_slice(&self.interp[slot * self.nnode..(slot + 1) * self.nnode]);

        if self.mixed {
            let dg_stride = self.geo_nnode * self.geo_ndim;
            deriv_geo
                .as_mut_data()
                .copy_from_slice(&self.deriv_geo[slot * dg_stride..(slot + 1) * dg_stride]);
            let d_stride = self.nnode * self.geo_ndim;
            deriv
                .as_mut_data()
                .copy_from_slice(&self.deriv[slot * d_stride..(slot + 1) * d_stride]);
        } else {
            let d_stride = self.nnode * self.geo_ndim;
            deriv
                .as_mut_data()
                .copy_from_slice(&self.deriv[slot * d_stride..(slot + 1) * d_stride]);
        }

        let jac_stride = self.space_ndim * self.geo_ndim;
        jacobian
            .as_mut_data()
            .copy_from_slice(&self.jacobian[slot * jac_stride..(slot + 1) * jac_stride]);

        let inv_stride = self.space_ndim * self.space_ndim;
        inv_jacobian
            .as_mut_data()
            .copy_from_slice(&self.inv_jacobian[slot * inv_stride..(slot + 1) * inv_stride]);

        let grad_stride = self.nnode * self.space_ndim;
        gradient
            .as_mut_data()
            .copy_from_slice(&self.gradient[slot * grad_stride..(slot + 1) * grad_stride]);

        self.det_jac[slot]
    }

    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    #[cfg(test)]
    pub(crate) fn cursor(&self) -> usize {
        self.cursor
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::MAX_CACHE_POINTS;
    use crate::integ::Gauss;
    use crate::shapes::scratchpad_testing::aux;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::mat_mat_mul;

    fn assert_send<T: Send>() {}

    fn assert_slice_bits_eq(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len());
        for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            assert_eq!(x.to_bits(), y.to_bits(), "bit mismatch at index {i}: {x} vs {y}");
        }
    }

    /// Evaluates `pad` at `ksi` bypassing the cache entirely to build a bit-exact reference
    fn uncached_gradient(reference: &mut Scratchpad, ksi: &[f64]) -> f64 {
        reference.calc_interp(ksi);
        let det_jac = reference.compute_jacobian_uncached(ksi).unwrap();
        if reference.geo_kind != reference.kind {
            (reference.fn_deriv)(&mut reference.deriv, ksi);
        }
        mat_mat_mul(
            &mut reference.gradient,
            1.0,
            &reference.deriv,
            &reference.inv_jacobian,
            0.0,
        )
        .unwrap();
        det_jac
    }

    fn set_distorted_qua8(pad: &mut Scratchpad, dx_mid: f64) {
        let coords = [
            [0.2, 0.1],
            [2.4, 0.3],
            [2.6, 1.9],
            [0.1, 1.7],
            [1.3 + dx_mid, 0.15],
            [2.55, 1.05],
            [1.35, 1.85],
            [0.12, 0.92],
        ];
        for (m, [x, y]) in coords.iter().enumerate() {
            pad.set_xx(m, 0, *x);
            pad.set_xx(m, 1, *y);
        }
    }

    #[test]
    fn scratchpad_is_send() {
        assert_send::<Scratchpad>();
    }

    #[test]
    fn multi_point_hits_are_bit_identical_for_single_and_mixed_pads() {
        let gauss = Gauss::new(GeoKind::Qua8);

        // Single-order distorted Qua8 across multiple passes and out-of-order lookups
        let mut pad = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        set_distorted_qua8(&mut pad, 0.05);
        let mut reference = pad.clone();
        assert!(pad.cache.is_none());

        for pass in 0..4 {
            for p in 0..gauss.npoint() {
                let iota = gauss.coords(p);
                let det_got = pad.calc_gradient(iota).unwrap();
                let det_ref = uncached_gradient(&mut reference, iota);
                assert_eq!(det_got.to_bits(), det_ref.to_bits(), "pass {pass}, point {p}");
                assert_slice_bits_eq(pad.interp.as_data(), reference.interp.as_data());
                assert_slice_bits_eq(pad.deriv.as_data(), reference.deriv.as_data());
                assert_slice_bits_eq(pad.jacobian.as_data(), reference.jacobian.as_data());
                assert_slice_bits_eq(pad.inv_jacobian.as_data(), reference.inv_jacobian.as_data());
                assert_slice_bits_eq(pad.gradient.as_data(), reference.gradient.as_data());
            }
            let cache = pad.cache.as_ref().unwrap();
            assert_eq!(cache.len(), gauss.npoint());
            assert_eq!(cache.cursor(), 0);
        }

        // Out-of-order query also hits without adding an extra slot
        for p in [4, 0, 8, 2, 7] {
            let iota = gauss.coords(p);
            let det_got = pad.calc_gradient(iota).unwrap();
            let det_ref = uncached_gradient(&mut reference, iota);
            assert_eq!(det_got.to_bits(), det_ref.to_bits());
            assert_slice_bits_eq(pad.gradient.as_data(), reference.gradient.as_data());
        }
        assert_eq!(pad.cache.as_ref().unwrap().len(), gauss.npoint());

        // Mixed-order Qua4 field over distorted Qua8 geometry
        let mut mixed = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua8).unwrap();
        set_distorted_qua8(&mut mixed, -0.07);
        let mut mixed_ref = mixed.clone();
        for _ in 0..3 {
            for p in 0..gauss.npoint() {
                let iota = gauss.coords(p);
                let det_got = mixed.calc_gradient(iota).unwrap();
                let det_ref = uncached_gradient(&mut mixed_ref, iota);
                assert_eq!(det_got.to_bits(), det_ref.to_bits());
                assert_slice_bits_eq(mixed.interp.as_data(), mixed_ref.interp.as_data());
                assert_slice_bits_eq(mixed.deriv.as_data(), mixed_ref.deriv.as_data());
                assert_slice_bits_eq(mixed.deriv_geo.as_data(), mixed_ref.deriv_geo.as_data());
                assert_slice_bits_eq(mixed.jacobian.as_data(), mixed_ref.jacobian.as_data());
                assert_slice_bits_eq(mixed.inv_jacobian.as_data(), mixed_ref.inv_jacobian.as_data());
                assert_slice_bits_eq(mixed.gradient.as_data(), mixed_ref.gradient.as_data());
            }
        }
        assert_eq!(mixed.cache.as_ref().unwrap().len(), gauss.npoint());
    }

    #[test]
    fn jacobian_hit_upgrades_to_gradient_without_clobbering_interp() {
        let mut mixed = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua8).unwrap();
        set_distorted_qua8(&mut mixed, 0.04);
        let mut reference = mixed.clone();
        let ksi = [0.25, -0.5];

        // Seed cache via calc_jacobian, then set recognizable sentinels in interp and gradient
        let dj1 = mixed.calc_jacobian(&ksi).unwrap();
        assert!(!mixed.cache.as_ref().unwrap().has_gradient(0));
        mixed.interp.fill(123.0);
        mixed.gradient.fill(456.0);
        mixed.deriv.fill(789.0);

        // Second calc_jacobian hit must NOT clobber interp, gradient, or mixed deriv
        let dj2 = mixed.calc_jacobian(&ksi).unwrap();
        assert_eq!(dj1.to_bits(), dj2.to_bits());
        assert_eq!(mixed.interp[0].to_bits(), 123.0_f64.to_bits());
        assert_eq!(mixed.gradient.get(0, 0).to_bits(), 456.0_f64.to_bits());
        assert_eq!(mixed.deriv.get(0, 0).to_bits(), 789.0_f64.to_bits());

        // Upgrading via calc_gradient fills interp, mixed deriv, and gradient bit-for-bit
        let dj3 = mixed.calc_gradient(&ksi).unwrap();
        let dj_ref = uncached_gradient(&mut reference, &ksi);
        assert_eq!(dj3.to_bits(), dj_ref.to_bits());
        assert!(mixed.cache.as_ref().unwrap().has_gradient(0));
        assert_slice_bits_eq(mixed.interp.as_data(), reference.interp.as_data());
        assert_slice_bits_eq(mixed.deriv.as_data(), reference.deriv.as_data());
        assert_slice_bits_eq(mixed.deriv_geo.as_data(), reference.deriv_geo.as_data());
        assert_slice_bits_eq(mixed.gradient.as_data(), reference.gradient.as_data());

        // Subsequent calc_gradient hits preserve all values
        mixed.gradient.fill(-999.0);
        let dj4 = mixed.calc_gradient(&ksi).unwrap();
        assert_eq!(dj4.to_bits(), dj_ref.to_bits());
        assert_slice_bits_eq(mixed.gradient.as_data(), reference.gradient.as_data());
    }

    #[test]
    fn set_xx_invalidates_cache_and_clone_resets_cache() {
        let mut pad = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        set_distorted_qua8(&mut pad, 0.0);
        let ksi = [0.0, 0.0];
        let dj_before = pad.calc_gradient(&ksi).unwrap();
        assert_eq!(pad.cache.as_ref().unwrap().len(), 1);

        // Clone starts with an empty cache and works independently
        let mut cloned = pad.clone();
        assert!(cloned.cache.is_none());
        let dj_clone = cloned.calc_gradient(&ksi).unwrap();
        assert_eq!(dj_clone.to_bits(), dj_before.to_bits());

        // Modifying node 4 via set_xx invalidates the cache immediately
        pad.set_xx(4, 0, 1.8);
        assert_eq!(pad.ok_xxt, false);
        assert_eq!(pad.cache.as_ref().unwrap().len(), 0);
        assert!(pad.calc_jacobian(&ksi).is_err());

        // Restoring ok_xxt via the last coordinate recomputes with the updated geometry
        pad.set_xx(7, 1, 1.35);
        assert_eq!(pad.ok_xxt, true);
        let dj_after = pad.calc_gradient(&ksi).unwrap();
        let mut reference = pad.clone();
        let dj_ref = uncached_gradient(&mut reference, &ksi);
        assert_eq!(dj_after.to_bits(), dj_ref.to_bits());
        assert_slice_bits_eq(pad.gradient.as_data(), reference.gradient.as_data());
        assert_ne!(dj_before.to_bits(), dj_after.to_bits());

        // First pad's clone stays unaffected
        let dj_clone_still = cloned.calc_gradient(&ksi).unwrap();
        assert_eq!(dj_clone_still.to_bits(), dj_before.to_bits());
    }

    #[test]
    fn degenerate_error_is_not_cached_and_cache_is_bounded_at_64() {
        // Degenerate zero-area triangle must fail repeatedly without caching
        let mut deg = Scratchpad::new(2, GeoKind::Tri3).unwrap();
        for m in 0..3 {
            deg.set_xx(m, 0, 1.0);
            deg.set_xx(m, 1, 2.0);
        }
        assert!(deg.calc_gradient(&[0.2, 0.3]).is_err());
        assert!(deg.calc_gradient(&[0.2, 0.3]).is_err());
        assert!(deg.cache.is_none());

        // Cable in 3D via calc_jacobian uses the cache without inv_jacobian/gradient
        let mut cable = aux::gen_scratchpad_with_coords(3, GeoKind::Lin3);
        let n1 = cable.calc_jacobian(&[-0.5]).unwrap();
        let n2 = cable.calc_jacobian(&[-0.5]).unwrap();
        assert_eq!(n1.to_bits(), n2.to_bits());

        // Evaluating > 64 points stays bounded at MAX_CACHE_POINTS
        let mut pad = aux::gen_scratchpad_with_coords(2, GeoKind::Qua8);
        let mut reference = pad.clone();
        for i in 0..80 {
            let t = -0.9 + 1.8 * (i as f64) / 79.0;
            let ksi = [t, 0.25];
            let got = pad.calc_gradient(&ksi).unwrap();
            let expected = uncached_gradient(&mut reference, &ksi);
            assert_eq!(got.to_bits(), expected.to_bits());
            assert_slice_bits_eq(pad.gradient.as_data(), reference.gradient.as_data());
        }
        assert_eq!(pad.cache.as_ref().unwrap().len(), MAX_CACHE_POINTS);

        // Re-probing one of the kept recent points still hits and stays bit-exact
        let recent = [-0.9 + 1.8 * (70.0 / 79.0), 0.25];
        let got = pad.calc_gradient(&recent).unwrap();
        let expected = uncached_gradient(&mut reference, &recent);
        assert_eq!(got.to_bits(), expected.to_bits());
        assert_slice_bits_eq(pad.gradient.as_data(), reference.gradient.as_data());
    }
}

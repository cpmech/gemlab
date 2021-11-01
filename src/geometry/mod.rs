#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circle {
    pub center: [f64; 2],
    pub radius: f64,
    pub tolerance: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointAndHalo {
    pub ndim: usize,      // 2D or 3D
    pub x: [[f64; 3]; 9], // 5=4+1 in 2D or 9=8+1 in 3D (each contains 3 coords)
}

pub(crate) fn compute_halo(halo: &mut [[f64; 3]; 9], ndim: usize, x: &[f64], dx: &[f64]) {
    assert!(ndim == 2 || ndim == 3);
    assert_eq!(x.len(), ndim);
    assert_eq!(dx.len(), ndim);
    if ndim == 2 {
        halo[0][0] = x[0] - dx[0];
        halo[0][1] = x[1] - dx[1];

        halo[1][0] = x[0] + dx[0];
        halo[1][1] = x[1] - dx[1];

        halo[2][0] = x[0] + dx[0];
        halo[2][1] = x[1] + dx[1];

        halo[3][0] = x[0] - dx[0];
        halo[3][1] = x[1] + dx[1];
    } else {
        halo[0][0] = x[0];
        halo[0][1] = x[1];
        halo[0][2] = x[2];

        halo[1][0] = x[0] - dx[0];
        halo[1][1] = x[1] - dx[1];
        halo[1][2] = x[2] - dx[2];

        halo[2][0] = x[0] + dx[0];
        halo[2][1] = x[1] - dx[1];
        halo[2][2] = x[2] - dx[2];

        halo[3][0] = x[0] + dx[0];
        halo[3][1] = x[1] + dx[1];
        halo[3][2] = x[2] - dx[2];

        halo[4][0] = x[0] - dx[0];
        halo[4][1] = x[1] + dx[1];
        halo[4][2] = x[2] - dx[2];

        halo[5][0] = x[0] - dx[0];
        halo[5][1] = x[1] - dx[1];
        halo[5][2] = x[2] + dx[2];

        halo[6][0] = x[0] + dx[0];
        halo[6][1] = x[1] - dx[1];
        halo[6][2] = x[2] + dx[2];

        halo[7][0] = x[0] + dx[0];
        halo[7][1] = x[1] + dx[1];
        halo[7][2] = x[2] + dx[2];

        halo[8][0] = x[0] - dx[0];
        halo[8][1] = x[1] + dx[1];
        halo[8][2] = x[2] + dx[2];
    }
}

impl PointAndHalo {
    pub fn new(x: &[f64], dx: &[f64]) -> Self {
        let ndim = x.len();
        assert_eq!(dx.len(), ndim);
        assert!(ndim >= 2 && ndim <= 3);
        let mut p = PointAndHalo { ndim, x: [[0.0; 3]; 9] };
        compute_halo(&mut p.x, ndim, x, dx);
        p
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circle_traits_work() {
        let circle = Circle {
            center: [-1.0, -1.0],
            radius: 2.0,
            tolerance: 1e-3,
        };
        let clone = circle.clone();
        assert_eq!(circle, clone);
        assert_eq!(
            format!("{:?}", circle),
            "Circle { center: [-1.0, -1.0], radius: 2.0, tolerance: 0.001 }"
        );
    }

    #[test]
    fn point_and_halo_traits_work() {
        let p = PointAndHalo {
            ndim: 2,
            x: [[0.0; 3]; 9],
        };
        let clone = p.clone();
        assert_eq!(p, clone);
        assert_eq!(
            format!("{:?}", p),
            "PointAndHalo { ndim: 2, x: [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] }"
        );
    }
}

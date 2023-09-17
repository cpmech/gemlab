use super::{Mesh, PointId};
use crate::shapes::Scratchpad;

impl Mesh {
    /// Sets the pad's matrix of coordinates X given a list of point ids
    ///
    /// # Panics
    ///
    /// 1. Make sure `pad.kind.nnode() == points.len()`; otherwise a panic will occur
    /// 2. This function does not check for bounds on point indices and dimensions
    /// 3. Use [Mesh::check_all] to capture (some) errors
    #[inline]
    pub fn set_pad(&self, pad: &mut Scratchpad, points: &[PointId]) {
        let nnode = pad.kind.nnode();
        assert_eq!(nnode, points.len());
        for m in 0..nnode {
            for j in 0..self.ndim {
                pad.set_xx(m, j, self.points[points[m]].coords[j]);
            }
        }
    }

    /// Returns the (min,max) point coordinates in a mesh
    #[inline]
    pub fn get_limits(&self) -> (Vec<f64>, Vec<f64>) {
        let mut min = vec![f64::MAX; self.ndim];
        let mut max = vec![f64::MIN; self.ndim];
        for point in &self.points {
            for i in 0..self.ndim {
                min[i] = f64::min(min[i], point.coords[i]);
                max[i] = f64::max(max[i], point.coords[i]);
            }
        }
        (min, max)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::mesh::Samples;
    use crate::shapes::Scratchpad;

    #[test]
    fn set_pad_coords_works_2d() {
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        let mesh = Samples::two_qua4();
        let cell = &mesh.cells[0];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad, &cell.points);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌         ┐\n\
             │ 0 1 1 0 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        let cell = &mesh.cells[1];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad, &cell.points);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌         ┐\n\
             │ 1 2 2 1 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
    }

    #[test]
    fn set_pad_coords_works_3d() {
        //       8-------------11
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  9-------------10    |
        //  |    .         |    |
        //  |    4---------|----7
        //  |   /.         |   /|
        //  |  / .         |  / |
        //  | /  .         | /  |
        //  |/   .         |/   |
        //  5--------------6    |
        //  |    .         |    |
        //  |    0---------|----3
        //  |   /          |   /
        //  |  /           |  /
        //  | /            | /
        //  |/             |/
        //  1--------------2
        let mesh = Samples::two_hex8();
        let cell = &mesh.cells[0];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad, &cell.points);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 0 0 0 0 1 1 1 1 │\n\
             └                 ┘"
        );
        let cell = &mesh.cells[1];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad, &cell.points);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 1 1 1 1 2 2 2 2 │\n\
             └                 ┘"
        );
    }

    #[test]
    fn get_mesh_limits_works() {
        let mesh = &Samples::two_qua4();
        let (min, max) = mesh.get_limits();
        assert_eq!(min, &[0.0, 0.0]);
        assert_eq!(max, &[2.0, 1.0]);

        let mesh = &Samples::two_hex8();
        let (min, max) = mesh.get_limits();
        assert_eq!(min, &[0.0, 0.0, 0.0]);
        assert_eq!(max, &[1.0, 1.0, 2.0]);
    }
}

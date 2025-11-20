use gemlab::prelude::*;
use gemlab::StrError;

fn main() -> Result<(), StrError> {
    let mesh = Mesh::new_zero_homogeneous(2, 4, 2, GeoKind::Tri3).unwrap();
    assert_eq!(mesh.ndim, 2);
    assert_eq!(mesh.points.len(), 4);
    assert_eq!(mesh.cells.len(), 2);
    for i in 0..4 {
        let point = &mesh.points[i];
        assert_eq!(point.id, i);
        assert_eq!(point.coords, vec![0.0, 0.0]);
    }
    for e in 0..2 {
        let cell = &mesh.cells[e];
        assert_eq!(cell.id, e);
        assert_eq!(cell.marker, 1);
        assert_eq!(cell.kind, GeoKind::Tri3);
        assert_eq!(cell.points, vec![0, 0, 0]);
    }
    assert_eq!(
        mesh.check_all().err(),
        Some("cannot compute inverse due to zero determinant")
    );
    Ok(())
}

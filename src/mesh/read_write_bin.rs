//! Binary mesh files carry a four-byte `GML1` magic and format-version header followed by a
//! `postcard` payload of the whole [Mesh].
//!
//! The layout is native 64-bit, and structural invariants are validated when reading.

use super::Mesh;
use crate::StrError;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

impl Mesh {
    /// Magic number and format version of binary mesh files
    const MAGIC: &'static [u8; 4] = b"GML1";

    /// Writes a binary mesh file
    ///
    /// The file starts with the four-byte `GML1` magic and format-version header followed by a
    /// `postcard` payload of the whole mesh.
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::{Mesh, Samples};
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     let mesh = Samples::two_qua4();
    ///     let path = "/tmp/gemlab/doc_mesh_bin.bin";
    ///     mesh.write_bin(&path)?;
    ///     let mesh_read = Mesh::read_bin(&path)?;
    ///     assert_eq!(mesh_read.points.len(), 6);
    ///     assert_eq!(mesh_read.cells.len(), 2);
    ///     assert_eq!(mesh_read.points[0].coords, vec![0.0, 0.0]);
    ///     Ok(())
    /// }
    /// ```
    pub fn write_bin<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        if self.cells.is_empty() {
            return Err("there are no cells to write");
        }
        let path = Path::new(full_path);
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let file = File::create(path).map_err(|_| "cannot create file")?;
        let mut writer = BufWriter::new(file);
        writer.write_all(Self::MAGIC).map_err(|_| "cannot write file")?;
        postcard::to_io(self, &mut writer).map_err(|_| "cannot write file")?;
        writer.flush().map_err(|_| "cannot flush file")?;
        Ok(())
    }

    /// Reads a binary mesh file
    ///
    /// The file must start with the four-byte `GML1` magic and format-version header followed by
    /// a `postcard` payload of the whole mesh; structural invariants are validated after decoding.
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::{Mesh, Samples};
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     let mesh = Samples::two_qua4();
    ///     let path = "/tmp/gemlab/doc_mesh_bin_read.bin";
    ///     mesh.write_bin(&path)?;
    ///     let mesh_read = Mesh::read_bin(&path)?;
    ///     assert_eq!(mesh_read.points.len(), 6);
    ///     assert_eq!(mesh_read.cells.len(), 2);
    ///     assert_eq!(mesh_read.points[0].coords, vec![0.0, 0.0]);
    ///     Ok(())
    /// }
    /// ```
    pub fn read_bin<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path);
        let bytes = fs::read(path).map_err(|_| "cannot open file")?;
        if !bytes.starts_with(Self::MAGIC) {
            return Err("cannot parse binary mesh file");
        }
        let (mesh, rest) = postcard::take_from_bytes::<Self>(&bytes[Self::MAGIC.len()..])
            .map_err(|_| "cannot parse binary mesh file")?;
        if !rest.is_empty() {
            return Err("cannot parse binary mesh file");
        }
        check_structure(&mesh)?;
        Ok(mesh)
    }
}

/// Checks the structural invariants of a mesh decoded from a binary file
fn check_structure(mesh: &Mesh) -> Result<(), StrError> {
    if mesh.ndim != 2 && mesh.ndim != 3 {
        return Err("ndim must be 2 or 3");
    }
    let npoint = mesh.points.len();
    if npoint < 2 {
        return Err("npoint must be ≥ 2");
    }
    if mesh.cells.is_empty() {
        return Err("ncell must be ≥ 1");
    }
    for point in &mesh.points {
        if point.coords.len() != mesh.ndim {
            return Err("point coords length must equal ndim");
        }
    }
    for cell in &mesh.cells {
        if cell.kind.ndim() > mesh.ndim {
            return Err("cell geo_ndim cannot be greater than ndim");
        }
    }
    mesh.check_ids_and_kind()?;
    for (_, p1, p2) in &mesh.marked_edges {
        if *p1 >= npoint || *p2 >= npoint {
            return Err("marked edge point id is out of bounds");
        }
    }
    for (_, p1, p2, p3, p4) in &mesh.marked_faces {
        if *p1 >= npoint || *p2 >= npoint || *p3 >= npoint || (*p4 != usize::MAX && *p4 >= npoint) {
            return Err("marked face point id is out of bounds");
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Mesh;
    use crate::mesh::{Cell, GeoKind, Point, Samples};
    use std::fs;

    fn check_same(a: &Mesh, b: &Mesh) {
        assert_eq!(a.ndim, b.ndim);
        assert_eq!(a.points.len(), b.points.len());
        for (pa, pb) in a.points.iter().zip(&b.points) {
            assert_eq!(pa.id, pb.id);
            assert_eq!(pa.marker, pb.marker);
            assert_eq!(pa.coords.len(), pb.coords.len());
            for i in 0..pa.coords.len() {
                assert_eq!(pa.coords[i].to_bits(), pb.coords[i].to_bits());
            }
        }
        assert_eq!(a.cells.len(), b.cells.len());
        for (ca, cb) in a.cells.iter().zip(&b.cells) {
            assert_eq!(ca.id, cb.id);
            assert_eq!(ca.marker, cb.marker);
            assert_eq!(ca.kind, cb.kind);
            assert_eq!(ca.points, cb.points);
        }
        assert_eq!(a.marked_edges, b.marked_edges);
        assert_eq!(a.marked_faces, b.marked_faces);
    }

    fn write_raw(mesh: &Mesh, path: &str) {
        let payload = postcard::to_stdvec(mesh).unwrap();
        let mut bytes = Vec::with_capacity(4 + payload.len());
        bytes.extend_from_slice(b"GML1");
        bytes.extend_from_slice(&payload);
        if let Some(p) = std::path::Path::new(path).parent() {
            fs::create_dir_all(p).unwrap();
        }
        fs::write(path, bytes).unwrap();
    }

    fn mesh_2d_two_points() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point {
                    id: 0,
                    marker: 0,
                    coords: vec![0.0, 0.0],
                },
                Point {
                    id: 1,
                    marker: 0,
                    coords: vec![1.0, 0.0],
                },
            ],
            cells: vec![Cell {
                id: 0,
                marker: 1,
                kind: GeoKind::Lin2,
                points: vec![0, 1],
            }],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        }
    }

    fn mesh_2d_four_points() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point {
                    id: 0,
                    marker: 0,
                    coords: vec![0.0, 0.0],
                },
                Point {
                    id: 1,
                    marker: 0,
                    coords: vec![1.0, 0.0],
                },
                Point {
                    id: 2,
                    marker: 0,
                    coords: vec![1.0, 1.0],
                },
                Point {
                    id: 3,
                    marker: 0,
                    coords: vec![0.0, 1.0],
                },
            ],
            cells: vec![Cell {
                id: 0,
                marker: 1,
                kind: GeoKind::Qua4,
                points: vec![0, 1, 2, 3],
            }],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        }
    }

    #[test]
    fn write_and_read_two_qua4() {
        let path = "/tmp/gemlab/test_bin_two_qua4.bin";
        let mesh = Samples::two_qua4();
        mesh.write_bin(path).unwrap();
        let mesh_read = Mesh::read_bin(path).unwrap();
        check_same(&mesh, &mesh_read);
    }

    #[test]
    fn write_and_read_two_hex8() {
        let path = "/tmp/gemlab/test_bin_two_hex8.bin";
        let mesh = Samples::two_hex8();
        mesh.write_bin(path).unwrap();
        let mesh_read = Mesh::read_bin(path).unwrap();
        check_same(&mesh, &mesh_read);
    }

    #[test]
    fn write_and_read_mixed_shapes_2d() {
        let path = "/tmp/gemlab/test_bin_mixed_shapes_2d.bin";
        let mut mesh = Mesh::read("./data/meshes/mixed_shapes_2d.msh").unwrap();
        mesh.marked_edges = vec![(-100, 0, 1), (-200, 1, 2)];
        mesh.write_bin(path).unwrap();
        let mesh_read = Mesh::read_bin(path).unwrap();
        check_same(&mesh, &mesh_read);
    }

    #[test]
    fn write_and_read_mixed_shapes_3d() {
        let path = "/tmp/gemlab/test_bin_mixed_shapes_3d.bin";
        let mut mesh = Mesh::read("./data/meshes/mixed_shapes_3d.msh").unwrap();
        mesh.marked_edges = vec![(-100, 0, 1)];
        mesh.marked_faces = vec![(-200, 0, 1, 2, usize::MAX), (-300, 0, 1, 2, 3)];
        mesh.write_bin(path).unwrap();
        let mesh_read = Mesh::read_bin(path).unwrap();
        check_same(&mesh, &mesh_read);
    }

    #[test]
    fn write_bin_captures_errors() {
        let mesh_empty = Mesh {
            ndim: 2,
            points: Samples::two_qua4().points,
            cells: Vec::new(),
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        assert_eq!(
            mesh_empty.write_bin("/tmp/gemlab/test_bin_no_cells.bin").err(),
            Some("there are no cells to write")
        );

        let mesh = Samples::two_qua4();
        assert_eq!(mesh.write_bin("/tmp/").err(), Some("cannot create file"));
        assert_eq!(
            mesh.write_bin("/dev/null/dir/test.bin").err(),
            Some("cannot create directory")
        );
        assert_eq!(
            Mesh::read_bin("/tmp/gemlab/test_bin_missing.bin").err(),
            Some("cannot open file")
        );
    }

    #[test]
    fn read_bin_rejects_bad_container() {
        fs::create_dir_all("/tmp/gemlab").unwrap();
        let mesh = Samples::two_qua4();
        let payload = postcard::to_stdvec(&mesh).unwrap();

        let path_bad_magic = "/tmp/gemlab/test_bin_bad_magic.bin";
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"XXXX");
        bytes.extend_from_slice(&payload);
        fs::write(path_bad_magic, bytes).unwrap();
        assert_eq!(
            Mesh::read_bin(path_bad_magic).err(),
            Some("cannot parse binary mesh file")
        );

        let path_short_magic = "/tmp/gemlab/test_bin_short_magic.bin";
        fs::write(path_short_magic, b"GML").unwrap();
        assert_eq!(
            Mesh::read_bin(path_short_magic).err(),
            Some("cannot parse binary mesh file")
        );

        let path_valid = "/tmp/gemlab/test_bin_container_valid.bin";
        mesh.write_bin(path_valid).unwrap();
        let valid = fs::read(path_valid).unwrap();

        let path_truncated = "/tmp/gemlab/test_bin_truncated.bin";
        fs::write(path_truncated, &valid[..valid.len() - 1]).unwrap();
        assert_eq!(
            Mesh::read_bin(path_truncated).err(),
            Some("cannot parse binary mesh file")
        );

        let path_extra = "/tmp/gemlab/test_bin_extra_byte.bin";
        let mut extra = valid.clone();
        extra.push(0);
        fs::write(path_extra, extra).unwrap();
        assert_eq!(Mesh::read_bin(path_extra).err(), Some("cannot parse binary mesh file"));
    }

    #[test]
    fn read_bin_rejects_bad_ndim() {
        for ndim in [1, 4] {
            let mesh = Mesh {
                ndim,
                points: vec![
                    Point {
                        id: 0,
                        marker: 0,
                        coords: vec![0.0; ndim],
                    },
                    Point {
                        id: 1,
                        marker: 0,
                        coords: vec![1.0; ndim],
                    },
                ],
                cells: vec![Cell {
                    id: 0,
                    marker: 1,
                    kind: GeoKind::Lin2,
                    points: vec![0, 1],
                }],
                marked_edges: Vec::new(),
                marked_faces: Vec::new(),
            };
            let path = format!("/tmp/gemlab/test_bin_bad_ndim_{}.bin", ndim);
            write_raw(&mesh, &path);
            assert_eq!(Mesh::read_bin(&path).err(), Some("ndim must be 2 or 3"));
        }
    }

    #[test]
    fn read_bin_rejects_bad_npoint() {
        let mesh = Mesh {
            ndim: 2,
            points: vec![Point {
                id: 0,
                marker: 0,
                coords: vec![0.0, 0.0],
            }],
            cells: vec![Cell {
                id: 0,
                marker: 1,
                kind: GeoKind::Lin2,
                points: vec![0, 0],
            }],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let path = "/tmp/gemlab/test_bin_bad_npoint.bin";
        write_raw(&mesh, path);
        assert_eq!(Mesh::read_bin(path).err(), Some("npoint must be ≥ 2"));
    }

    #[test]
    fn read_bin_rejects_bad_ncell() {
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point {
                    id: 0,
                    marker: 0,
                    coords: vec![0.0, 0.0],
                },
                Point {
                    id: 1,
                    marker: 0,
                    coords: vec![1.0, 0.0],
                },
            ],
            cells: Vec::new(),
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let path = "/tmp/gemlab/test_bin_bad_ncell.bin";
        write_raw(&mesh, path);
        assert_eq!(Mesh::read_bin(path).err(), Some("ncell must be ≥ 1"));
    }

    #[test]
    fn read_bin_rejects_bad_coords() {
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point {
                    id: 0,
                    marker: 0,
                    coords: vec![0.0, 0.0],
                },
                Point {
                    id: 1,
                    marker: 0,
                    coords: vec![1.0, 0.0, 0.0],
                },
            ],
            cells: vec![Cell {
                id: 0,
                marker: 1,
                kind: GeoKind::Lin2,
                points: vec![0, 1],
            }],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let path = "/tmp/gemlab/test_bin_bad_coords.bin";
        write_raw(&mesh, path);
        assert_eq!(Mesh::read_bin(path).err(), Some("point coords length must equal ndim"));
    }

    #[test]
    fn read_bin_rejects_bad_geo_ndim() {
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point {
                    id: 0,
                    marker: 0,
                    coords: vec![0.0, 0.0],
                },
                Point {
                    id: 1,
                    marker: 0,
                    coords: vec![1.0, 0.0],
                },
                Point {
                    id: 2,
                    marker: 0,
                    coords: vec![1.0, 1.0],
                },
                Point {
                    id: 3,
                    marker: 0,
                    coords: vec![0.0, 1.0],
                },
            ],
            cells: vec![Cell {
                id: 0,
                marker: 1,
                kind: GeoKind::Tet4,
                points: vec![0, 1, 2, 3],
            }],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let path = "/tmp/gemlab/test_bin_bad_geo_ndim.bin";
        write_raw(&mesh, path);
        assert_eq!(
            Mesh::read_bin(path).err(),
            Some("cell geo_ndim cannot be greater than ndim")
        );
    }

    #[test]
    fn read_bin_rejects_bad_ids() {
        let path = "/tmp/gemlab/test_bin_bad_point_id.bin";
        let mut mesh = mesh_2d_two_points();
        mesh.points[1].id = 2;
        write_raw(&mesh, path);
        assert_eq!(
            Mesh::read_bin(path).err(),
            Some("incorrect point id found; ids must be sequential")
        );

        let path = "/tmp/gemlab/test_bin_bad_cell_id.bin";
        let mut mesh = mesh_2d_two_points();
        mesh.cells[0].id = 2;
        write_raw(&mesh, path);
        assert_eq!(
            Mesh::read_bin(path).err(),
            Some("incorrect cell id found; ids must be sequential")
        );

        let path = "/tmp/gemlab/test_bin_bad_cell_nnode.bin";
        let mut mesh = mesh_2d_two_points();
        mesh.cells[0].kind = GeoKind::Qua4;
        write_raw(&mesh, path);
        assert_eq!(
            Mesh::read_bin(path).err(),
            Some("number of cell points does not correspond to kind.nnode()")
        );

        let path = "/tmp/gemlab/test_bin_bad_cell_point.bin";
        let mut mesh = mesh_2d_two_points();
        mesh.cells[0].points = vec![0, 2];
        write_raw(&mesh, path);
        assert_eq!(
            Mesh::read_bin(path).err(),
            Some("the id of a point specified in the points list of a cell is out of bounds")
        );
    }

    #[test]
    fn read_bin_rejects_bad_marked_edges() {
        let mut mesh = mesh_2d_two_points();
        mesh.marked_edges.push((-1, 0, 2));
        let path = "/tmp/gemlab/test_bin_bad_marked_edge.bin";
        write_raw(&mesh, path);
        assert_eq!(
            Mesh::read_bin(path).err(),
            Some("marked edge point id is out of bounds")
        );
    }

    #[test]
    fn read_bin_rejects_bad_marked_faces() {
        let bad_faces = [
            (1, 4, 1, 2, usize::MAX),
            (1, 0, 4, 2, usize::MAX),
            (1, 0, 1, 4, usize::MAX),
            (1, 0, 1, 2, 4),
        ];
        for (k, face) in bad_faces.iter().enumerate() {
            let mut mesh = mesh_2d_four_points();
            mesh.marked_faces.push(*face);
            let path = format!("/tmp/gemlab/test_bin_bad_marked_face_{}.bin", k);
            write_raw(&mesh, &path);
            assert_eq!(
                Mesh::read_bin(&path).err(),
                Some("marked face point id is out of bounds")
            );
        }

        let mut mesh = mesh_2d_four_points();
        mesh.marked_faces.push((1, 0, 1, 2, usize::MAX));
        let path = "/tmp/gemlab/test_bin_good_marked_face.bin";
        write_raw(&mesh, path);
        assert!(Mesh::read_bin(path).is_ok());
    }
}

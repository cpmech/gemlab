use super::Mesh;
use crate::StrError;
use std::ffi::OsStr;
use std::fmt::Write;
use std::fs::{self, File};
use std::io::Write as IoWrite;
use std::path::Path;

impl Mesh {
    /// Writes Paraview's VTU file
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write_vtu<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let ndim = self.ndim;
        let npoint = self.points.len();
        let ncell = self.cells.len();
        if ncell < 1 {
            return Err("there are no cells to write");
        }

        let mut buffer = String::new();

        // header
        write!(
            &mut buffer,
            "<?xml version=\"1.0\"?>\n\
         <VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\
         <UnstructuredGrid>\n\
         <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n",
            npoint, ncell
        )
        .unwrap();

        // nodes: coordinates
        write!(
            &mut buffer,
            "<Points>\n\
         <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n",
        )
        .unwrap();
        for index in 0..npoint {
            for dim in 0..ndim {
                write!(&mut buffer, "{:?} ", self.points[index].coords[dim]).unwrap();
            }
            if ndim == 2 {
                write!(&mut buffer, "0.0 ").unwrap();
            }
        }
        write!(
            &mut buffer,
            "\n</DataArray>\n\
         </Points>\n"
        )
        .unwrap();

        // elements: connectivity
        write!(
            &mut buffer,
            "<Cells>\n\
         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n"
        )
        .unwrap();
        for cell in &self.cells {
            if cell.kind.vtk_type().is_none() {
                return Err("cannot generate VTU file because VTK cell type is not available");
            }
            for p in &cell.points {
                write!(&mut buffer, "{} ", p).unwrap();
            }
        }

        // elements: offsets
        write!(
            &mut buffer,
            "\n</DataArray>\n\
         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
        )
        .unwrap();
        let mut offset = 0;
        for cell in &self.cells {
            offset += cell.points.len();
            write!(&mut buffer, "{} ", offset).unwrap();
        }

        // elements: types
        write!(
            &mut buffer,
            "\n</DataArray>\n\
         <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
        )
        .unwrap();
        for cell in &self.cells {
            if let Some(vtk) = cell.kind.vtk_type() {
                write!(&mut buffer, "{} ", vtk).unwrap();
            }
        }
        write!(
            &mut buffer,
            "\n</DataArray>\n\
         </Cells>\n"
        )
        .unwrap();

        write!(
            &mut buffer,
            "</Piece>\n\
         </UnstructuredGrid>\n\
         </VTKFile>\n"
        )
        .unwrap();

        // create directory
        let path = Path::new(full_path);
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }

        // write file
        let mut file = File::create(path).map_err(|_| "cannot create file")?;
        file.write_all(buffer.as_bytes()).map_err(|_| "cannot write file")?;

        // force sync
        file.sync_all().map_err(|_| "cannot sync file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{mesh::Samples, StrError};
    use std::fs;

    #[test]
    fn write_vtu_handles_unavailable_types() {
        assert_eq!(
            Samples::lin_cells().write_vtu("/tmp/gemlab/nothing.vtu").err(),
            Some("cannot generate VTU file because VTK cell type is not available")
        );
        assert_eq!(
            Samples::tri_cells().write_vtu("/tmp/gemlab/nothing.vtu").err(),
            Some("cannot generate VTU file because VTK cell type is not available")
        );
        assert_eq!(
            Samples::qua_cells().write_vtu("/tmp/gemlab/nothing.vtu").err(),
            Some("cannot generate VTU file because VTK cell type is not available")
        );
    }

    #[test]
    fn write_vtu_works_qua8_tri6_lin2() -> Result<(), StrError> {
        let mesh = Samples::qua8_tri6_lin2();
        let file_path = "/tmp/gemlab/test_qua8_tri6_lin2.vtu";
        mesh.write_vtu(file_path)?;
        let contents = fs::read_to_string(file_path).map_err(|_| "cannot open file")?;
        assert_eq!(
            contents,
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
<UnstructuredGrid>
<Piece NumberOfPoints="11" NumberOfCells="4">
<Points>
<DataArray type="Float64" NumberOfComponents="3" format="ascii">
0.0 0.0 0.0 0.5 0.0 0.0 1.0 0.0 0.0 1.433 0.25 0.0 1.866 0.5 0.0 1.433 0.75 0.0 1.0 1.0 0.0 0.5 1.0 0.0 0.0 1.0 0.0 0.0 0.5 0.0 1.0 0.5 0.0 
</DataArray>
</Points>
<Cells>
<DataArray type="Int32" Name="connectivity" format="ascii">
0 2 6 8 1 10 7 9 2 4 6 3 5 10 2 10 10 6 
</DataArray>
<DataArray type="Int32" Name="offsets" format="ascii">
8 14 16 18 
</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">
23 22 3 3 
</DataArray>
</Cells>
</Piece>
</UnstructuredGrid>
</VTKFile>
"#
        );
        Ok(())
    }

    #[test]
    fn write_vtu_works_mixed_shapes_3d() -> Result<(), StrError> {
        let mesh = Samples::mixed_shapes_3d();
        let file_path = "/tmp/gemlab/test_mixed_shapes_3d.vtu";
        mesh.write_vtu(file_path)?;
        let contents = fs::read_to_string(file_path).map_err(|_| "cannot open file")?;
        assert_eq!(
            contents,
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
<UnstructuredGrid>
<Piece NumberOfPoints="13" NumberOfCells="5">
<Points>
<DataArray type="Float64" NumberOfComponents="3" format="ascii">
0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1.0 0.0 1.0 1.0 1.0 1.0 0.0 1.0 1.0 1.0 2.0 0.0 0.0 2.0 0.0 0.0 2.0 1.0 1.0 -0.5 0.0 1.0 -1.0 0.0 
</DataArray>
</Points>
<Cells>
<DataArray type="Int32" Name="connectivity" format="ascii">
0 1 2 3 4 5 6 7 2 8 3 6 3 9 10 7 8 9 3 1 12 11 
</DataArray>
<DataArray type="Int32" Name="offsets" format="ascii">
8 12 16 19 22 
</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">
12 10 9 5 21 
</DataArray>
</Cells>
</Piece>
</UnstructuredGrid>
</VTKFile>
"#
        );
        Ok(())
    }
}

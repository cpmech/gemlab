use gemlab::mesh::Mesh;
use gemlab::mesh::Unstructured;
use gemlab::StrError;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::path::PathBuf;
use structopt::StructOpt;

const OUT_DIR: &str = "/tmp/gemlab";

/// Command line options
#[derive(Debug, StructOpt)]
#[structopt(name = "msh2tri", about = "Generates triangular meshes from a mesh file.")]
struct Options {
    /// Input file
    #[structopt(parse(from_os_str))]
    input: PathBuf,

    /// Generate Tri6 instead of Tri3 elements
    #[structopt(short, long)]
    tri6: bool,

    /// Apply node renumbering to reduce the matrix bandwidth
    #[structopt(short, long)]
    renumber: bool,

    /// Write a binary mesh file instead of the ASCII MSH file
    #[structopt(long)]
    bin: bool,

    #[structopt(long)]
    global_max_area: Option<f64>,
}

/// Reads a mesh file in either the ASCII MSH or the binary format
fn read_mesh(in_path: &Path) -> Result<Mesh, StrError> {
    let mut file = File::open(in_path).map_err(|_| "cannot open file")?;
    let mut magic = [0u8; 4];
    if file.read_exact(&mut magic).is_ok() && magic == *b"GML1" {
        Mesh::read_bin(in_path)
    } else {
        Mesh::read(in_path)
    }
}

fn main() -> Result<(), StrError> {
    // parse options
    let options = Options::from_args();

    // extract file stem
    let in_path = Path::new(&options.input);
    let fn_stem = in_path
        .file_stem()
        .ok_or("cannot get file stem")?
        .to_str()
        .ok_or("cannot convert file stem to str")?;

    // read mesh file as PSLG
    let pslg = read_mesh(in_path)?;

    // data for Trigen
    let holes = Vec::new();
    let o2 = options.tri6;
    let max_areas = None;
    let global_min_angle = None;
    let renumber = options.renumber;

    // call Trigen
    let mesh = Unstructured::call_trigen(
        &pslg,
        &holes,
        o2,
        max_areas,
        options.global_max_area,
        global_min_angle,
        renumber,
    )?;

    // save mesh file
    if options.bin {
        let path_out = format!("{}/{}.bin", OUT_DIR, fn_stem);
        mesh.write_bin(&path_out)?;
        println!("Generated binary mesh file: {}", path_out);
    } else {
        let path_out = format!("{}/{}.msh", OUT_DIR, fn_stem);
        mesh.write(&path_out)?;
        println!("Generated MSH file: {}", path_out);
    }
    Ok(())
}

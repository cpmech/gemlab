use gemlab::mesh::Mesh;
use gemlab::mesh::Unstructured;
use gemlab::StrError;
use std::path::Path;
use std::path::PathBuf;
use structopt::StructOpt;

const OUT_DIR: &str = "/tmp/gemlab";

/// Command line options
#[derive(Debug, StructOpt)]
#[structopt(name = "msh2tri", about = "Generates triangular meshes from a MSH file.")]
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

    #[structopt(long)]
    global_max_area: Option<f64>,
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

    // read MSH file as PSLG
    let pslg = Mesh::read(in_path)?;

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

    // save MSH file
    let path_out = format!("{}/{}.msh", OUT_DIR, fn_stem);
    mesh.write(&path_out)?;
    println!("Generated MSH file: {}", path_out);
    Ok(())
}

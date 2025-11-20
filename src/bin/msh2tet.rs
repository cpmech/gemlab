use gemlab::mesh::Mesh;
use gemlab::mesh::Unstructured;
use gemlab::StrError;
use std::path::Path;
use std::path::PathBuf;
use structopt::StructOpt;

// TODO: Finish msh2tet implementation

const OUT_DIR: &str = "/tmp/gemlab";

/// Command line options
#[derive(Debug, StructOpt)]
#[structopt(name = "msh2tet", about = "Generates tetrahedral meshes from a MSH file.")]
struct Options {
    /// Input file
    #[structopt(parse(from_os_str))]
    input: PathBuf,

    /// Point defining the region
    #[structopt(
        long,
        value_names = &["x", "y", "z"],
        number_of_values = 3,
        allow_hyphen_values = true
    )]
    region: Vec<f64>,

    /// Generate Tet10 instead of Tet4 elements
    #[structopt(short, long)]
    tet10: bool,

    /// Apply node renumbering to reduce the matrix bandwidth
    #[structopt(short, long)]
    renumber: bool,

    #[structopt(long)]
    global_max_volume: Option<f64>,
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

    // read MSH file as PLC
    let plc = Mesh::read(in_path)?;

    // data for Trigen
    let regions = vec![(1, options.region[0], options.region[1], options.region[2])];
    let holes = Vec::new();
    let o2 = options.tet10;
    let max_volumes = None;
    let global_min_angle = None;
    let renumber = options.renumber;

    // call Trigen
    let mesh = Unstructured::call_tetgen(
        &plc,
        &regions,
        &holes,
        o2,
        max_volumes,
        options.global_max_volume,
        global_min_angle,
        renumber,
    )?;

    // save MSH file
    let path_out = format!("{}/{}.msh", OUT_DIR, fn_stem);
    mesh.write(&path_out)?;
    println!("Generated MSH file: {}", path_out);
    Ok(())
}

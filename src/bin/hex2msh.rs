use gemlab::mesh::{Blocks3d, Draw, GeoKind, Structured};
use gemlab::StrError;
use std::path::Path;
use std::path::PathBuf;
use structopt::StructOpt;

/// Command line options
#[derive(Debug, StructOpt)]
#[structopt(
    name = "hex2msh",
    about = "Generates hex meshes from on a set of blocks. Outputs MSH and VTU files (VTU: only for n8 and n20)."
)]
struct Options {
    /// Input JSON file
    #[structopt(parse(from_os_str))]
    input: PathBuf,

    /// Output directory
    out_dir: String,

    /// Number of nodes on the target element (8, 20, 32)
    #[structopt(short, long, default_value = "8")]
    nnode: usize,

    /// Generate SVG figure with the wireframe of the mesh
    #[structopt(short, long)]
    svg_figure: bool,

    /// Show points in the SVG figure
    #[structopt(short = "p", long)]
    show_points: bool,

    /// Apply node renumbering to reduce the matrix bandwidth
    #[structopt(short, long)]
    renumber: bool,
}

fn main() -> Result<(), StrError> {
    // parse options
    let options = Options::from_args();

    // load input data from JSON file
    let in_path = Path::new(&options.input);
    let blocks = Blocks3d::read_json(&in_path)?;
    let fn_stem = in_path
        .file_stem()
        .ok_or("cannot get file stem")?
        .to_str()
        .ok_or("cannot convert file stem to str")?;

    // generate mesh
    let (target, do_vtu) = match options.nnode {
        8 => (GeoKind::Hex8, true),
        20 => (GeoKind::Hex20, true),
        32 => (GeoKind::Hex32, false),
        _ => return Err("invalid number of nodes for hexahedral cell"),
    };
    let mesh = Structured::from_blocks_3d(&blocks, target, options.renumber)?;

    // write MSH file
    let path_msh = format!("{}/{}.msh", options.out_dir, fn_stem);
    mesh.write(&path_msh)?;
    println!("\nGenerated MSH file: {}/{}.msh", options.out_dir, fn_stem);

    // write VTU file
    if do_vtu {
        let path_vtu = format!("{}/{}.vtu", options.out_dir, fn_stem);
        mesh.write_vtu(&path_vtu)?;
        println!("Generated VTU file: {}/{}.vtu", options.out_dir, fn_stem);
    }

    // write SVG file with wireframe
    if options.svg_figure {
        let mut draw = Draw::new();
        draw.size(1000.0, 1000.0);
        draw.show_point_dots(options.show_points);
        draw.all(&mesh, &format!("{}/{}.svg", options.out_dir, fn_stem))?;
        println!("Generated SVG file: {}/{}.svg", options.out_dir, fn_stem);
    }
    println!();
    Ok(())
}

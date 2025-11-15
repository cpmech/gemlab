use gemlab::mesh::{Draw, Mesh};
use gemlab::StrError;
use std::path::Path;
use std::path::PathBuf;
use structopt::StructOpt;

const OUT_DIR: &str = "/tmp/gemlab";

/// Command line options
#[derive(Debug, StructOpt)]
#[structopt(name = "drawmsh", about = "Draw MSH file")]
struct Options {
    /// Input MSH file
    #[structopt(parse(from_os_str))]
    input: PathBuf,

    /// Show points in the SVG figure
    #[structopt(short, long)]
    dots: bool,

    /// Show IDs of points in the SVG figure
    #[structopt(short, long)]
    point_ids: bool,

    /// Show IDs of cells in the SVG figure
    #[structopt(short, long)]
    cell_ids: bool,

    #[structopt(short, long)]
    width: Option<f64>,

    #[structopt(short, long)]
    height: Option<f64>,
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

    // load MSH file
    let mesh = Mesh::read(in_path)?;

    // write SVG file
    let mut fig = Draw::new();
    fig.show_point_dots(options.dots)
        .show_point_ids(options.point_ids)
        .show_cell_ids(options.cell_ids)
        .size(options.width.unwrap_or(800.0), options.height.unwrap_or(800.0));
    fig.all(&mesh, &format!("{}/{}.svg", OUT_DIR, fn_stem))?;
    println!("Generated SVG file: {}/{}.svg", OUT_DIR, fn_stem);
    Ok(())
}

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

    /// Skip mesh checks
    #[structopt(short = "k", long)]
    skip_check: bool,

    /// Show points in the SVG figure
    #[structopt(short, long)]
    dots: bool,

    /// Show IDs of points in the SVG figure
    #[structopt(short, long)]
    point_ids: bool,

    /// Show IDs of cells in the SVG figure
    #[structopt(short, long)]
    cell_ids: bool,

    /// Show point markers
    #[structopt(short = "v", long)]
    point_markers: bool,

    /// Show cell attributes
    #[structopt(short = "a", long)]
    cell_att: bool,

    /// Show edge markers
    #[structopt(short, long)]
    edge_markers: bool,

    /// Show face markers
    #[structopt(short, long)]
    face_markers: bool,

    /// Show normal vectors on boundaries
    #[structopt(short, long)]
    normals: bool,

    /// Hide cells
    #[structopt(short = "i", long)]
    hide_cells: bool,

    /// Boundary edges
    #[structopt(short = "y", long)]
    boundary_edges: bool,

    /// Boundary faces
    #[structopt(short = "z", long)]
    boundary_faces: bool,

    /// Figure width in points
    #[structopt(short, long)]
    width: Option<f64>,

    /// Figure height in points
    #[structopt(short, long)]
    height: Option<f64>,

    /// View the figure interactively in a window (also saves SVG file)
    #[structopt(long)]
    view: bool,

    /// Multiplier for drawing area range
    #[structopt(long)]
    m_range: Option<f64>,

    /// Multiplier for normal vector length
    #[structopt(long)]
    m_normal_vector: Option<f64>,

    /// Multiplier for normal vector marker length
    #[structopt(long)]
    m_normal_vector_marker: Option<f64>,

    /// Camera elevation angle in degrees (3D only)
    #[structopt(long)]
    elevation: Option<f64>,

    /// Camera azimuth angle in degrees (3D only)
    #[structopt(long)]
    azimuth: Option<f64>,

    /// Enable or disable glyph indicating XYZ directions (3D only)
    #[structopt(long = "glyph-enabled")]
    glyph_enabled: Option<bool>,

    /// Glyph origin (x y z) (3D only)
    #[structopt(
        long,
        value_names = &["x", "y", "z"],
        number_of_values = 3,
        allow_hyphen_values = true
    )]
    glyph_origin: Option<Vec<f64>>,

    /// Glyph size (3D only)
    #[structopt(long)]
    glyph_size: Option<f64>,
}

fn get_glyph_origin(values: &Option<Vec<f64>>) -> Result<[f64; 3], StrError> {
    match values {
        Some(v) if v.len() == 3 => Ok([v[0], v[1], v[2]]),
        Some(_) => Err("glyph-origin expects exactly three values"),
        None => Ok([0.0, 0.0, 0.0]),
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

    // load MSH file
    let mesh = Mesh::read(in_path)?;

    // check mesh data
    if !options.skip_check {
        mesh.check_all()?;
    }

    // configure drawing
    let mut draw = Draw::new();

    // Set boolean flags
    draw.show_cells(!options.hide_cells)
        .show_point_dots(options.dots)
        .show_point_ids(options.point_ids)
        .show_cell_ids(options.cell_ids)
        .show_point_marker(options.point_markers)
        .show_cell_att(options.cell_att)
        .show_edge_markers(options.edge_markers)
        .show_face_markers(options.face_markers)
        .show_normal_vectors(options.normals)
        .show_boundary_edges_3d(options.boundary_edges)
        .show_boundary_faces(options.boundary_faces)
        .set_view_flag(options.view);

    // Set figure size
    if options.width.is_some() || options.height.is_some() {
        draw.set_size(options.width.unwrap_or(800.0), options.height.unwrap_or(800.0));
    }

    // Set multipliers
    if let Some(m_range) = options.m_range {
        draw.set_m_range(m_range);
    }
    if let Some(m_normal_vector) = options.m_normal_vector {
        draw.set_m_normal_vector(m_normal_vector);
    }
    if let Some(m_normal_vector_marker) = options.m_normal_vector_marker {
        draw.set_m_normal_vector_marker(m_normal_vector_marker);
    }

    // Set camera for 3D meshes
    if mesh.ndim == 3 {
        let elevation = options.elevation.unwrap_or(30.0);
        let azimuth = options.azimuth.unwrap_or(30.0);
        draw.set_camera(elevation, azimuth);

        if let Some(enabled) = options.glyph_enabled {
            draw.show_glyph_3d(enabled);
        }

        if options.glyph_origin.is_some() || options.glyph_size.is_some() {
            let origin = get_glyph_origin(&options.glyph_origin)?;
            let size = options.glyph_size.unwrap_or(1.0);
            draw.set_glyph_3d(origin[0], origin[1], origin[2], size);
        }
    }

    // Write SVG file
    draw.all(&mesh, &format!("{}/{}.svg", OUT_DIR, fn_stem))?;
    println!("Generated SVG file: {}/{}.svg", OUT_DIR, fn_stem);
    Ok(())
}

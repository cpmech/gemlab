use gemlab::mesh::{draw_mesh, Mesh};
use gemlab::StrError;

fn main() -> Result<(), StrError> {
    // 1.0  5------,6.------7
    //      | [3],'   `.[4] |
    //      |  ,'       `.  |
    //      |,'           `.|
    // 0.5  3      [2]      4
    //      |`.           .'|
    //      |  `.       .'  |
    //      | [0]`.   .'[1] |
    // 0.0  0------`1'------2
    //     0.0     0.5     1.0
    let mesh = Mesh::from_text(
        r"# header
          # ndim npoint ncell
               2      8     5
          
          # points
          # id marker x y
             0 0 0.0 0.0
             1 0 0.5 0.0
             2 0 1.0 0.0
             3 0 0.0 0.5
             4 0 1.0 0.5
             5 0 0.0 1.0
             6 0 0.5 1.0
             7 0 1.0 1.0
          
          # cells
          # id attribute kind point_ids...
             0 1 tri3 0 1 3
             1 1 tri3 1 2 4
             2 1 qua4 1 4 6 3
             3 1 tri3 5 3 6
             4 1 tri3 7 6 4
         ",
    )?;
    println!("{}", mesh);
    assert_eq!(mesh.points.len(), 8);
    assert_eq!(mesh.cells.len(), 5);
    assert_eq!(mesh.cells[0].points.len(), 3);
    assert_eq!(mesh.cells[2].points.len(), 4);
    draw_mesh(&mesh, true, false, false, "/tmp/gemlab/example_mesh_2d_tri3_qua4.svg")?;
    Ok(())
}

use gemlab::shapes::GeoKind;

fn main() {
    let kind = GeoKind::Qua8;
    for t in 0..kind.triangulate_ntriangle() {
        for i in 0..3 {
            let m = kind.triangulate_triangle_nodes(t, i);
            let ksi = if m >= kind.nnode() {
                let k = m - kind.nnode();
                kind.triangulate_extra_coords(k)
            } else {
                kind.reference_coords(m)
            };
            println!("t = {t}, i = {i}, m = {m}, ksi = ({}, {})", ksi[0], ksi[1]);
        }
    }
}

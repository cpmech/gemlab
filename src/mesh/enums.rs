/// Defines the location of points
pub enum At {
    /// Fixed x, any (y,z). 2D (vertical-line) or 3D (yz-plane)
    X(f64),

    /// Fixed y, any (x,z). 2D (horizontal-line) or 3D (xz-plane)
    Y(f64),

    /// Fixed z, any (y,z). 3D only (xy-plane)
    Z(f64),

    /// Fixed (x,y), any z. 2D (point) or 3D (z-line)
    XY(f64, f64),

    /// Fixed (y,z), any x. 3D only (x-line)
    YZ(f64, f64),

    /// Fixed (x,z), any y. 3D only (y-line)
    XZ(f64, f64),

    /// At (x,y,z). 3D only (point)
    XYZ(f64, f64, f64),

    /// Circumference of circle at (x,y,radius). 2D only
    Circle(f64, f64, f64),

    /// Surface of cylinder given two axes. 3D only
    ///
    /// Holds: (axis_1_x, axis_1_y, axis_1_z, axis_2_x, axis_2_y, axis_2_z, radius).
    Cylinder(f64, f64, f64, f64, f64, f64, f64),
}

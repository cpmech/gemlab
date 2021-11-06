/// Defines the location of points
pub enum At {
    /// At (x,y)
    XY(f64, f64),

    /// At (x,y,z)
    XYZ(f64, f64, f64),

    /// Horizontal line passing at (x,y)
    Horizontal(f64, f64),

    /// Vertical line passing at (x,y)
    Vertical(f64, f64),

    /// Circumference of circle at (x,y,radius)
    Circle(f64, f64, f64),

    /// Circumference of circle at (x,y,z,radius)
    Circle3D(f64, f64, f64, f64),

    /// Surface of cylinder parallel to x (x,y,z,radius)
    CylinderX(f64, f64, f64, f64),

    /// Surface of cylinder parallel to y (x,y,z,radius)
    CylinderY(f64, f64, f64, f64),

    /// Surface of cylinder parallel to z (x,y,z,radius)
    CylinderZ(f64, f64, f64, f64),

    /// At plane with normal along x and passing at (x,y,z)
    PlaneNormalX(f64, f64, f64),

    /// At plane with normal along y and passing at (x,y,z)
    PlaneNormalY(f64, f64, f64),

    /// At plane with normal along z and passing at (x,y,z)
    PlaneNormalZ(f64, f64, f64),
}

/// Defines geometric features by locating individual points or points along edges or faces
pub enum Geo {
    Point(At),
    Edge(At),
    Face(At),
}

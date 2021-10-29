# Idea

```text
#[derive(Clone, Copy, Debug)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Vector2D {
    pub u: f64,
    pub v: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Vector3D {
    pub u: f64,
    pub v: f64,
    pub w: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Circle {
    pub center: Point2D,
    pub radius: f64,
    pub tolerance: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Cylinder {
    pub axis: Vector3D,
    pub center: Point3D,
    pub radius: f64,
    pub tolerance: f64,
}
```

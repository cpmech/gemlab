use crate::shapes::GeoKind;

/// Defines the interface to access Cell/Face information
pub trait AsCell {
    /// Returns the geometry kind
    fn kind(&self) -> GeoKind;

    /// Returns the marker
    fn marker(&self) -> i32;

    /// Returns an access to the list of point ids
    fn points(&self) -> &[usize];
}

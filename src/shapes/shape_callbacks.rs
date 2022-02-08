use super::{i32_to_fn_deriv, i32_to_fn_interp, GeoKind};
use russell_lab::{Matrix, Vector};
use serde::{de::Deserialize, de::Deserializer, Serializer};
use std::fmt;

/// Defines an alias for interpolation functions
#[derive(Clone)]
pub(crate) struct FnInterp(pub(crate) GeoKind, pub(crate) fn(&mut Vector, &[f64]));

/// Defines an alias for derivative of interpolation functions
#[derive(Clone)]
pub(crate) struct FnDeriv(pub(crate) GeoKind, pub(crate) fn(&mut Matrix, &[f64]));

/// Implements the serialize method for FnInterp
/// Only the GeoKind as i32 needs to be encoded
pub(crate) fn shape_serialize_fn_interp<S>(interp: &FnInterp, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_i32(interp.0 as i32)
}

/// Implements the deserialize method for FnInterp
/// Decodes the GeoKind as i32 and find the proper function
pub(crate) fn shape_deserialize_fn_interp<'de, D>(deserializer: D) -> Result<FnInterp, D::Error>
where
    D: Deserializer<'de>,
{
    let kind: i32 = Deserialize::deserialize(deserializer)?;
    Ok(i32_to_fn_interp(kind))
}

/// Implements the serialize method for FnDeriv
/// Only the GeoKind as i32 needs to be encoded
pub(crate) fn shape_serialize_fn_deriv<S>(deriv: &FnDeriv, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_i32(deriv.0 as i32)
}

/// Implements the deserialize method for FnDeriv
/// Decodes the GeoKind as i32 and find the proper function
pub(crate) fn shape_deserialize_fn_deriv<'de, D>(deserializer: D) -> Result<FnDeriv, D::Error>
where
    D: Deserializer<'de>,
{
    let kind: i32 = Deserialize::deserialize(deserializer)?;
    Ok(i32_to_fn_deriv(kind))
}

impl fmt::Debug for FnInterp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "FnInterp {{ {} }}", self.0 as i32)
    }
}

impl fmt::Debug for FnDeriv {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "FnDeriv {{ {} }}", self.0 as i32)
    }
}

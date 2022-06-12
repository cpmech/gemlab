#![cfg_attr(coverage_nightly, feature(no_coverage))]

//! Gemlab -- Geometry, meshes, and integration for finite element analyses

/// Defines a type alias for the error type as a static string
pub type StrError = &'static str;

pub mod geometry;
pub mod integ;
pub mod mesh;
pub mod shapes;
pub mod util;

// run code from README file
#[cfg(doctest)]
mod test_readme {
    macro_rules! external_doc_test {
        ($x:expr) => {
            #[doc = $x]
            extern "C" {}
        };
    }
    external_doc_test!(include_str!("../README.md"));
}

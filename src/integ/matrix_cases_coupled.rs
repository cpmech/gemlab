use crate::StrError;

/// Implements the shape(N) time scalar(S) time gradient(G) integration case with different shapes (e.g., coupling matrix)
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠       →
/// Kᵐⁿ = │ NBᵐ s Gⁿ dΩ
///       ⌡
///       Ωₑ
/// ```
pub fn mat_coupling_nsg() -> Result<(), StrError> {
    Err("mat_coupling_nsg: TODO")
}

/// Implements the gradient(G) time tensor(T) time shape(N) integration case with different shapes (e.g., coupling matrix)
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠ →
/// Kᵐⁿ = │ GBᵐ ⋅ T Nⁿ dΩ
///       ⌡       ▔
///       Ωₑ
/// ```
pub fn mat_coupling_gtn() -> Result<(), StrError> {
    Err("mat_coupling_gtn: TODO")
}

/// Implements the shape(N) time vector(V) time shape(N) integration case with different shapes (e.g., coupling matrix)
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠    →
/// Kᵐⁿ = │ Nᵐ v NBⁿ dΩ
///       ⌡
///       Ωₑ
/// ```
pub fn mat_coupling_nvn() -> Result<(), StrError> {
    Err("mat_coupling_nvn: TODO")
}

/// Implements the gradient(G) time scalar(S) time shape(N) integration case with different shapes (e.g., coupling matrix)
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠ →
/// Kᵐⁿ = │ Gᵐ s NBⁿ dΩ
///       ⌡
///       Ωₑ
/// ```
pub fn mat_coupling_gsn() -> Result<(), StrError> {
    Err("mat_coupling_gsn: TODO")
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{mat_coupling_gsn, mat_coupling_gtn, mat_coupling_nsg, mat_coupling_nvn};

    #[test]
    fn functions_return_todo() {
        assert_eq!(mat_coupling_nsg().err(), Some("mat_coupling_nsg: TODO"));
        assert_eq!(mat_coupling_gtn().err(), Some("mat_coupling_gtn: TODO"));
        assert_eq!(mat_coupling_nvn().err(), Some("mat_coupling_nvn: TODO"));
        assert_eq!(mat_coupling_gsn().err(), Some("mat_coupling_gsn: TODO"));
    }
}

/// Defines the kind of shape
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Kind {
    Lin2,
    Lin3,
    Lin4,
    Lin5,
    Tri3,
    Tri6,
    Tri10,
    Tri15,
    Qua4,
    Qua8,
    Qua9,
    Qua12,
    Qua16,
    Tet4,
    Tet10,
    Hex8,
    Hex20,
}

/// Defines Qua shapes
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum KindQua {
    Qua4,
    Qua8,
    Qua9,
    Qua12,
    Qua16,
}

/// Defines Hex shapes
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum KindHex {
    Hex8,
    Hex20,
}

/// Defines Qua or Hex shapes
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum KindQuaOrHex {
    Qua4,
    Qua8,
    Qua9,
    Qua12,
    Qua16,
    Hex8,
    Hex20,
}

pub fn kind_from_ndim_npoint(ndim: usize, npoint: usize) -> Option<Kind> {
    match (ndim, npoint) {
        (1, 2) => Some(Kind::Lin2),
        (1, 3) => Some(Kind::Lin3),
        (1, 4) => Some(Kind::Lin4),
        (1, 5) => Some(Kind::Lin5),
        (2, 3) => Some(Kind::Tri3),
        (2, 6) => Some(Kind::Tri6),
        (2, 10) => Some(Kind::Tri10),
        (2, 15) => Some(Kind::Tri15),
        (2, 4) => Some(Kind::Qua4),
        (2, 8) => Some(Kind::Qua8),
        (2, 9) => Some(Kind::Qua9),
        (2, 12) => Some(Kind::Qua12),
        (2, 16) => Some(Kind::Qua16),
        (3, 4) => Some(Kind::Tet4),
        (3, 10) => Some(Kind::Tet10),
        (3, 8) => Some(Kind::Hex8),
        (3, 20) => Some(Kind::Hex20),
        _ => None,
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Kind, KindHex, KindQua};

    #[test]
    fn kind_enums_work() {
        let lin = Kind::Lin2;
        let clone = lin.clone();
        assert_eq!(lin, clone);
        assert_eq!(format!("{:?}", lin), "Lin2");

        let qua = KindQua::Qua12;
        let clone = qua.clone();
        assert_eq!(qua, clone);
        assert_eq!(format!("{:?}", qua), "Qua12");

        let hex = KindHex::Hex20;
        let clone = hex.clone();
        assert_eq!(hex, clone);
        assert_eq!(format!("{:?}", hex), "Hex20");
    }
}

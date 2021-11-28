use super::*;
use crate::StrError;

impl Shape {
    /// Selects a set of integrations points
    pub fn select_int_points(&mut self, nip: usize) -> Result<(), StrError> {
        match self.class {
            GeoClass::Lin => match nip {
                1 => self.ip_data = &IP_LIN_LEGENDRE_1,
                2 => self.ip_data = &IP_LIN_LEGENDRE_2,
                3 => self.ip_data = &IP_LIN_LEGENDRE_3,
                4 => self.ip_data = &IP_LIN_LEGENDRE_4,
                5 => self.ip_data = &IP_LIN_LEGENDRE_5,
                _ => return Err("integration points set is not available for this GeoClass"),
            },
            GeoClass::Tri => match nip {
                1 => self.ip_data = &IP_TRI_INTERNAL_1,
                3 => self.ip_data = &IP_TRI_INTERNAL_3,
                _ => return Err("integration points set is not available for this GeoClass"),
            },
            GeoClass::Qua => match nip {
                1 => self.ip_data = &IP_QUA_LEGENDRE_1,
                4 => self.ip_data = &IP_QUA_LEGENDRE_4,
                9 => self.ip_data = &IP_QUA_LEGENDRE_9,
                16 => self.ip_data = &IP_QUA_LEGENDRE_16,
                _ => return Err("integration points set is not available for this GeoClass"),
            },
            GeoClass::Tet => match nip {
                1 => self.ip_data = &IP_TET_INTERNAL_1,
                4 => self.ip_data = &IP_TET_INTERNAL_4,
                5 => self.ip_data = &IP_TET_INTERNAL_5,
                6 => self.ip_data = &IP_TET_INTERNAL_6,
                _ => return Err("integration points set is not available for this GeoClass"),
            },
            GeoClass::Hex => match nip {
                8 => self.ip_data = &IP_HEX_LEGENDRE_8,
                9 => self.ip_data = &IP_HEX_WILSON_CORNER_9,
                6 => self.ip_data = &IP_HEX_IRONS_6,
                14 => self.ip_data = &IP_HEX_IRONS_14,
                27 => self.ip_data = &IP_HEX_LEGENDRE_27,
                _ => return Err("integration points set is not available for this GeoClass"),
            },
        }
        Ok(())
    }

    pub fn integrate_fx(&mut self) -> Result<f64, StrError> {
        Ok(0.0)
    }
}

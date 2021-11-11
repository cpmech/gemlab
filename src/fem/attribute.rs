use super::ElementKind;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Attribute {
    pub(super) kind: ElementKind,
    pub(super) inactive: bool,
    pub(super) parameters: HashMap<String, f64>,
    pub(super) flags: HashMap<String, bool>,
}

impl Attribute {
    pub fn new(kind: ElementKind) -> Self {
        Attribute {
            kind,
            inactive: false,
            parameters: HashMap::new(),
            flags: HashMap::new(),
        }
    }

    pub fn set_parameter(&mut self, name: &str, value: f64) -> &mut Self {
        if let Some(parameter) = self.parameters.get_mut(name) {
            *parameter = value;
        }
        self
    }

    pub fn set_flag(&mut self, name: &str, value: bool) -> &mut Self {
        if let Some(flag) = self.flags.get_mut(name) {
            *flag = value;
        }
        self
    }
}

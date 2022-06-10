use gemlab::util::PI;
use gemlab::StrError;
use russell_lab::Vector;

struct Processor<'a, T> {
    callback: Option<Box<dyn Fn(&mut Vector, &[f64], &T) -> Result<(), StrError> + 'a>>,
}

impl<'a, T> Processor<'a, T> {
    fn set_callback(&mut self, c: impl Fn(&mut Vector, &[f64], &T) -> Result<(), StrError> + 'a) {
        self.callback = Some(Box::new(c));
    }

    fn process_events(&self, args: &T) -> Result<(), StrError> {
        let ksi = vec![0.0; 2];
        let mut x = Vector::new(2);
        match &self.callback {
            Some(cb) => cb(&mut x, &ksi, args)?,
            None => (),
        }
        println!("got x =\n{}", x);
        Ok(())
    }
}

fn simple_callback<T>(_: &mut Vector, _: &[f64], _: &T) -> Result<(), StrError> {
    Err("nothing to see here")
}

pub struct NaturalToReal {
    pub rmin: f64,
    pub rmax: f64,
    pub amin: f64,
    pub amax: f64,
    pub zmin: f64,
    pub zmax: f64,
}

impl NaturalToReal {
    pub fn new() -> Self {
        NaturalToReal {
            rmin: 5.0,
            rmax: 10.0,
            amin: 30.0 * PI / 180.0,
            amax: 60.0 * PI / 180.0,
            zmin: 0.0,
            zmax: 1.0,
        }
    }
}

pub fn natural_to_real_cylindrical(x: &mut Vector, ksi: &[f64], args: &NaturalToReal) -> Result<(), StrError> {
    if x.dim() != ksi.len() {
        return Err("x.dim() must be equal to ksi.len()");
    }
    const KSI_MIN: f64 = -1.0;
    const KSI_DEL: f64 = 2.0;
    let r = args.rmin + (ksi[0] - KSI_MIN) * (args.rmax - args.rmin) / KSI_DEL;
    let a = args.amin + (ksi[1] - KSI_MIN) * (args.amax - args.amin) / KSI_DEL;
    x[0] = r * f64::cos(a);
    x[1] = r * f64::sin(a);
    if x.dim() == 3 {
        x[2] = args.zmin + (ksi[2] - KSI_MIN) * (args.zmax - args.zmin) / KSI_DEL;
    }
    Ok(())
}

fn main() -> Result<(), StrError> {
    let mut p = Processor {
        callback: Some(Box::new(simple_callback)),
        // callback: None,
    };
    let args = NaturalToReal::new();
    println!("{:?}", p.process_events(&args).err());
    let v = Vector::from(&[1.0, 2.0]);
    let my_fun = move |x: &mut Vector, ksi: &[f64], _: &NaturalToReal| {
        println!("using ksi = {:?}", ksi);
        x[0] = v[0];
        x[1] = v[1];
        Ok(())
    };
    // println!("v = {:?}", v);
    p.set_callback(my_fun);
    p.process_events(&args)?;
    p.set_callback(natural_to_real_cylindrical);
    p.process_events(&args)?;
    Ok(())
}

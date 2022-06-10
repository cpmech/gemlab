use gemlab::StrError;
use russell_lab::Vector;

struct Processor<'a> {
    callback: Option<Box<dyn FnMut(&mut Vector, &[f64]) -> Result<(), StrError> + 'a>>,
}

impl<'a> Processor<'a> {
    fn set_callback(&mut self, c: impl FnMut(&mut Vector, &[f64]) -> Result<(), StrError> + 'a) {
        self.callback = Some(Box::new(c));
    }

    fn process_events(&mut self) -> Result<(), StrError> {
        let ksi = vec![0.0; 2];
        let mut x = Vector::new(2);
        // (self.callback)(&mut x, &ksi);
        match &mut self.callback {
            Some(a) => a(&mut x, &ksi)?,
            None => (),
        }
        println!("got x =\n{}", x);
        Ok(())
    }
}

fn simple_callback(_: &mut Vector, _: &[f64]) -> Result<(), StrError> {
    Err("nothing to see here")
}

fn main() -> Result<(), StrError> {
    let mut p = Processor {
        callback: Some(Box::new(simple_callback)),
    };
    println!("{:?}", p.process_events().err());
    let v = Vector::from(&[1.0, 2.0]);
    let my_fun = move |x: &mut Vector, ksi: &[f64]| {
        println!("using ksi = {:?}", ksi);
        x[0] = v[0];
        x[1] = v[1];
        Ok(())
    };
    // println!("v = {:?}", v);
    p.set_callback(my_fun);
    p.process_events()?;
    Ok(())
}

use russell_lab::Vector;

struct Processor<'a> {
    callback: Box<dyn FnMut(&mut Vector, &[f64]) + 'a>,
}

impl<'a> Processor<'a> {
    fn set_callback(&mut self, c: impl FnMut(&mut Vector, &[f64]) + 'a) {
        self.callback = Box::new(c);
    }

    fn process_events(&mut self) {
        let ksi = vec![0.0; 2];
        let mut x = Vector::new(2);
        (self.callback)(&mut x, &ksi);
        println!("got x =\n{}", x);
    }
}

fn simple_callback(x: &mut Vector, ksi: &[f64]) {
    println!("ksi = {:?}", ksi);
    println!("{}", x);
}

fn main() {
    let mut p = Processor {
        callback: Box::new(simple_callback),
    };
    p.process_events();
    let v = Vector::from(&[1.0, 2.0]);
    let x0 = 1.0;
    let x1 = 2.0;
    let my_fun = move |x: &mut Vector, ksi: &[f64]| {
        println!("using ksi = {:?}", ksi);
        x[0] = v[0];
        x[1] = v[1];
    };
    println!("x0 = {}", x0);
    println!("x1 = {}", x1);
    // println!("v = {:?}", v);
    p.set_callback(my_fun);
    p.process_events();
}

#[macro_use]
extern crate approx; // For the macro relative_eq!
extern crate nalgebra as na;
use na::{Vector3};

use mean_field_dipolar::get_dd_int;

fn main() {
    let dip_v = Vector3::new(0., 0., 1.);
    let dist_v = Vector3::new(1., 0., 0.);

    let interaction = get_dd_int(dist_v, dip_v);

    println!("d-d interaction {}", interaction);
}

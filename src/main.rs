extern crate nalgebra as na;
extern crate serde_json;
use mean_field_dipolar::dipolar::{simulation_sweep, Pattern};
use mean_field_dipolar::lattice::get_checkerboard;
use mean_field_dipolar::util;
use na::DVector;
use std::f64::consts::PI;
use std::path::Path;
use std::time::Instant;


fn main() {
    // save path
    let save_path = Path::new("./");

    let patt = Pattern::HStripe;

    simulation_sweep(save_path, &patt, (1, 2), (4, 6), PI/2., 0., 20.);
}

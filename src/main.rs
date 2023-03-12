extern crate nalgebra as na;
extern crate serde_json;
use mean_field_dipolar::dipolar::{simulation_sweep, Pattern};
use mean_field_dipolar::lattice::get_checkerboard;
use mean_field_dipolar::util;
use na::DVector;
use std::path::Path;


fn main() {
    // save path
    let save_path = Path::new("./");

    let patt = Pattern::CB;

    simulation_sweep(save_path, &patt, (1, 2), (4, 5), 0., 0., 20.);
}

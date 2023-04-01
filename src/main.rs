extern crate nalgebra as na;
extern crate serde_json;
extern crate toml;

use mean_field_dipolar::dipolar::{simulation_sweep, Pattern};
use mean_field_dipolar::lattice::get_checkerboard;
use mean_field_dipolar::util;
use na::DVector;
use std::f64::consts::PI;
use std::path::Path;
use std::time::Instant;
use std::fs;
use serde::Deserialize;


#[derive(Deserialize)]
struct Config {
    /// pattern to simulate
   pattern: String,
   /// interaction range start
   range_start: usize,
   /// interaction range end
   range_end: usize,
   /// size start
   size_start: usize,
   /// size end
   size_end: usize,
   /// theta (in fraction of PI)
   theta: f64,
   /// phi (in fraction of PI)
   phi: f64, 
   /// onsite interaction
   u_onsite: f64,
}

fn main() {
    // save path
    let save_path = Path::new("./");
    let config_path = Path::new("./sim.toml");

    // read the config
    let config_str = fs::read_to_string(config_path)
        .expect("Config file wasn't found");

    let config: Config = toml::from_str(config_str.as_ref())
        .expect("Failed reading config file");

    let patt = util::parse_pattern_str(config.pattern);

    simulation_sweep(save_path, &patt, 
                     (config.range_start, config.range_end+1),
                     (config.size_start, config.size_end+1),
                     config.theta*PI,
                     config.phi*PI, config.u_onsite);
}

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

    let patt = match config.pattern.to_lowercase().as_ref() {
        "filled" => Pattern::Filled,
        "cb" => Pattern::CB,
        "hstripe" => Pattern::HStripe,
        _ => Pattern::Filled
    };

    println!("Selected pattern {:?}", patt);

    simulation_sweep(save_path, &patt, (1, 2), (4, 6), PI/2., 0., 20.);
}

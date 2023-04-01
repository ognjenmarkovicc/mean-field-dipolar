extern crate nalgebra as na;
extern crate serde_json;
extern crate toml;

use mean_field_dipolar::dipolar::{simulation_sweep};
use mean_field_dipolar::util;
use std::f64::consts::PI;
use std::path::{Path, PathBuf};
use std::fs;
use clap::{Parser};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Simulation results path
    #[arg(short, long, value_name = "FILE")]
    respath: Option<PathBuf>,
}

fn main() {
    let cli = Cli::parse();

    // save path
    let save_path = match cli.respath.as_deref() {
        Some(path) => Path::new(path),
        None => Path::new("./")
    };
    println!("Selected results save path: {save_path:?}");

    // config file path
    let config_path = Path::new("./sim.toml");

    // read the config
    let config_str = fs::read_to_string(config_path)
        .expect("Config file wasn't found");

    let config: util::Config = toml::from_str(config_str.as_ref())
        .expect("Failed reading config file");

    let patt = util::parse_pattern_str(config.pattern);

    simulation_sweep(save_path, &patt, 
                     (config.range_start, config.range_end+1),
                     (config.size_start, config.size_end+1),
                     config.theta*PI,
                     config.phi*PI, config.u_onsite);
}

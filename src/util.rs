use na::{DVector, Scalar};
use std::fs;
use std::path::Path;
use serde::ser;
use serde::Deserialize;
use super::dipolar::Pattern;

/// Basic linspace function
/// 
/// Doesn't perform any checks, use with
/// care. Will silently convert num and idx to f64 
/// even if num cannot be represented as f64
pub fn linspace(start: f64, stop: f64, num: usize, endpoint: bool) -> DVector<f64> {

    let denom = if endpoint {
        num - 1
    } else {
        num
    };

    let delta = (stop - start)/(denom as f64);
    DVector::from_iterator(num, (0..num).map(|idx| (idx as f64)*delta + start))
}

/// Save a DVector<T> into a json file
pub fn save_vector_json<T, P>(filename: P, values: DVector<T>)
where P: AsRef<Path>, 
      T: Scalar + ser::Serialize,
{
    // Convert the vector to a JSON string
    let json_string = serde_json::to_string(&values).unwrap();                                         

    match fs::write(filename, json_string) {
        Ok(_) => (),
        Err(s) => println!("Error writing to file: {s}"),
    }   
}

/// Parse pattern from a config string
pub fn parse_pattern_str(pattern_str: String) -> Pattern {
    let patt = match pattern_str.to_lowercase().as_ref() {
        "filled" => Pattern::Filled,
        "cb" => Pattern::CB,
        "hstripe" => Pattern::HStripe,
        _ => Pattern::Filled
    };
    patt
}

/// Pattern struct
#[derive(Deserialize)]
pub struct Config {
   /// pattern to simulate
   pub pattern: String,
   /// interaction range start
   pub range_start: usize,
   /// interaction range end
   pub range_end: usize,
   /// size start
   pub size_start: usize,
   /// size end
   pub size_end: usize,
   /// theta (in fraction of PI)
   pub theta: f64,
   /// phi (in fraction of PI)
   pub phi: f64, 
   /// onsite interaction
   pub u_onsite: f64,
}
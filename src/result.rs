use na::{DVector, Scalar};
use std::fs;
use std::path::Path;
use serde::ser;

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

extern crate nalgebra as na;
extern crate serde_json;
use mean_field_dipolar::dipolar::{DipolarSystem, get_dd_int_site, generate_dd_int_mat, generate_mat_m, get_mu_inequality, get_tunneling};
use mean_field_dipolar::patterns::get_checkerboard;
use mean_field_dipolar::result;
use mean_field_dipolar::util::linspace;
use na::DVector;
use std::path::Path;


fn main() {
    // save path
    let save_path = Path::new("./");

    let mut dip_system = DipolarSystem::new(0., 0., 20., 1, 4);
    dip_system.update_occupation(get_checkerboard(&dip_system.latt));
    generate_dd_int_mat(&mut dip_system);
    let (lower, upper) = get_mu_inequality(&dip_system);

    println!("Lower mu {:.2} upper mu {:.2}", lower, upper);

    if lower < upper {
        let no_points = 100;
        let mu_vals = linspace(lower, upper, no_points, true);

        let tunneling = DVector::from_iterator(no_points, 
                                               mu_vals.iter()
                                                      .map(|mu| get_tunneling(*mu, &dip_system,
                                                                              1.,  1e-3, 1e-2)));

        result::save_vector_json(save_path.join("tunneling.json")
                                 , tunneling);
        result::save_vector_json(save_path.join("mu.json")
                                 , mu_vals);

    }

}

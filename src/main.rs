extern crate nalgebra as na;
use na::{Vector3};

use mean_field_dipolar::{get_dd_int, PeriodicLattice};

fn main() {
    let dip_v = Vector3::new(0., 0., 1.);
    let dist_v = Vector3::new(1., 0., 0.);

    let interaction = get_dd_int(dist_v, dip_v);

    println!("d-d interaction {}", interaction);

    let system = PeriodicLattice { system_size: 4 };
    let idx = 5;
    assert_eq!(1, system.get_idx_periodic(idx));
}

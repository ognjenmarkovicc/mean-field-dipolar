extern crate nalgebra as na;
use approx::relative_eq;
use na::{Vector3};

use mean_field_dipolar::lattice::{PeriodicLattice,
    SpinIdx, LattPos};
use mean_field_dipolar::dipolar::{DipolarSystem, get_dd_int, get_dd_int_site};
use mean_field_dipolar::patterns::get_checkerboard;

fn main() {
    let dip_v = Vector3::new(0., 0., 1.);
    let dist_v = Vector3::new(1., 0., 0.);

    let interaction = get_dd_int(dist_v, dip_v);

    println!("d-d interaction {}", interaction);

    let system = PeriodicLattice::new(4);
    let idx = 5;
    assert_eq!(1, system.get_idx_periodic(idx));

    let sp = SpinIdx::new(5, &system);

    let pos = LattPos::from(sp);

    println!("{:?}", pos);

    let mat = get_checkerboard(&system);

    println!("{:?}", mat);

    let dip_system = DipolarSystem {theta: 0., phi: 0.,
                                    u_onsite: 0., int_range: 2};
    let interaction = get_dd_int_site(0, 0, &dip_system,
                                      &mat, &system);

    println!("Interaction {}", interaction);

}

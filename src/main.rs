extern crate nalgebra as na;
use na::{Vector3};

use mean_field_dipolar::lattice::{PeriodicLattice,
    SpinIdx, LattPos};
use mean_field_dipolar::dipolar::get_dd_int;
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
}

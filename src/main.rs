extern crate nalgebra as na;
use approx::relative_eq;
use na::{Vector3};

use mean_field_dipolar::lattice::{PeriodicLattice,
    SpinIdx, LattPos};
use mean_field_dipolar::dipolar::{DipolarSystem, get_dd_int, get_dd_int_site, generate_dd_int_mat, generate_mat_m};
use mean_field_dipolar::patterns::get_checkerboard;

fn main() {
    let dip_v = Vector3::new(0., 0., 1.);
    let dist_v = Vector3::new(1., 0., 0.);

    let interaction = get_dd_int(dist_v, dip_v);

    println!("d-d interaction {}", interaction);

    let system = PeriodicLattice::new(2);
    let idx = 1;
    assert_eq!(1, system.get_idx_periodic(idx));

    let sp = SpinIdx::new(1, &system);

    let pos = LattPos::from(sp);

    println!("{:?}", pos);

    let mat = get_checkerboard(&system);

    println!("{:?}", mat);

    let dip_system = DipolarSystem::new(0., 0., 20., 1, 4);

    let interaction = get_dd_int_site(0, 0, &dip_system,
                                      &mat, &system);

    println!("Interaction {}", interaction);

    let dd_mat = generate_dd_int_mat(&dip_system, &mat, &system);

    println!("{:?}", dd_mat);

    let t = 1.;
    let mu = 1.;
    let m_mat = generate_mat_m(mu, t, &dip_system, &mat,
                               &system, &dd_mat);

    println!("{:?}", m_mat);
}

#[macro_use]
extern crate approx; // For the macro relative_eq!
extern crate nalgebra as na;
use na::Vector3;
use core::panic;
use std::f64::consts::PI;

pub mod lattice;
pub mod dipolar;
pub mod patterns;

#[cfg(test)]
mod tests {
    use crate::{patterns::get_checkerboard, dipolar::get_dd_int_site};

    use super::*;
    use dipolar::{get_dd_int, DipolarSystem};
    use lattice::{LattPos, PeriodicLattice, SpinIdx};

    #[test]
    fn dd_repulsive_test() {
        let dip_v = Vector3::new(0., 0., 1.);
        let dist_v = Vector3::new(1., 0., 0.);

        let int = get_dd_int(dist_v, dip_v);

        relative_eq!(int, 1.);
    }


    #[test]
    fn dd_attractive_test() {
        let dip_v = Vector3::new(1., 0., 0.);
        let dist_v = Vector3::new(1., 0., 0.);

        let int = get_dd_int(dist_v, dip_v);

        relative_eq!(int, -2.);
    }

    #[test]
    fn periodic_idx_pos_test() {
        let system = PeriodicLattice::new(4);
        let idx = 5;
        assert_eq!(1, system.get_idx_periodic(idx));
    }

    #[test]
    fn periodic_idx_neg_test() {
        let system = PeriodicLattice::new(4);
        let idx = -1;
        assert_eq!(3, system.get_idx_periodic(idx));
    }

    #[test]
    #[should_panic]
    fn latt_pos_panic_test() {
        let system = PeriodicLattice::new(4);
        let _ = LattPos::new(6, 5, &system);
    }

    #[test]
    fn pos_to_spin() {
        let system = PeriodicLattice::new(4);
        let pos = LattPos::new(3, 1, &system);
        let sp = SpinIdx::from(pos);

        assert_eq!(sp.idx, 7);
    }

    #[test]
    fn spin_to_pos() {
        let system = PeriodicLattice::new(4);
        let sp = SpinIdx::new(7, &system);
        let pos = LattPos::from(sp);

        assert_eq!(pos.x, 3);
        assert_eq!(pos.y, 1);
    }

    #[test]
    fn dipole_vec_test() {
        let dip_system = DipolarSystem {theta: PI/2., phi: 0., u_onsite: 0., int_range: 2};
        let dip_vec = dip_system.get_dipole_vec();

        relative_eq!(dip_vec[0], 1.);
    }

    #[test]
    fn dipole_int_test() {
        let dip_system = DipolarSystem {theta: 0., phi: 0., u_onsite: 0., int_range: 1};
        let system = PeriodicLattice::new(4);
        let occupation = get_checkerboard(&system);

        let interaction = get_dd_int_site(0, 0, &dip_system,
                                          &occupation, &system);
        assert_eq!(interaction, 4.);  
    }
}
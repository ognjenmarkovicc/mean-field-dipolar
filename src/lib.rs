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
    use crate::{patterns::get_checkerboard,
                dipolar::{get_dd_int_site, generate_mat_m, generate_dd_int_mat}};

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
        let dip_system = DipolarSystem::new(PI/2., 0., 0., 2, 4);

        let dip_vec = dip_system.get_dipole_vec();

        relative_eq!(dip_vec[0], 1.);
    }

    #[test]
    fn dipole_int_test() {
        let dip_system = DipolarSystem::new(0., 0., 0., 1, 4);

        let system = PeriodicLattice::new(4);
        let occupation = get_checkerboard(&system);

        let interaction = get_dd_int_site(0, 0, &dip_system,
                                          &occupation, &system);
        assert_eq!(interaction, 4.);  
    }

    #[test]
    fn m_matrix_det_test() {
        let dip_system = DipolarSystem::new(0., 0., 20., 1, 4);


        let system = PeriodicLattice::new(2);
        let occupation = get_checkerboard(&system);

        let dd_mat = generate_dd_int_mat(&dip_system, &occupation, &system);
        
        let t = 1.;
        let mu = 1.;
        let m_mat = generate_mat_m(mu, t, &dip_system, &occupation,
                                   &system, &dd_mat);
        
        // value taken from the python version of the
        assert_eq!(m_mat.determinant(), -0.4736842105263155);
    }
}
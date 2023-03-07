#[macro_use]
extern crate approx; // For the macro assert_relative_eq!
extern crate nalgebra as na;
use na::Vector3;
use std::f64::consts::PI;

pub mod lattice;
pub mod dipolar;
pub mod patterns;

#[cfg(test)]
mod tests {
    use crate::{patterns::get_checkerboard,
                dipolar::{get_dd_int_site, generate_mat_m, generate_dd_int_mat, get_tunneling, get_mu_inequality}};

    use super::*;
    use dipolar::{get_dd_int, DipolarSystem};
    use lattice::{LattPos, PeriodicLattice, SpinIdx};

    #[test]
    fn dd_repulsive_test() {
        let dip_v = Vector3::new(0., 0., 1.);
        let dist_v = Vector3::new(1., 0., 0.);

        let int = get_dd_int(dist_v, dip_v);

        assert_relative_eq!(int, 1.);
    }


    #[test]
    fn dd_attractive_test() {
        let dip_v = Vector3::new(1., 0., 0.);
        let dist_v = Vector3::new(1., 0., 0.);

        let int = get_dd_int(dist_v, dip_v);

        assert_relative_eq!(int, -2.);
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

        assert_relative_eq!(dip_vec[0], 1.);
    }

    #[test]
    fn dipole_int_test() {
        let mut dip_system = DipolarSystem::new(0., 0., 0., 1, 4);

        dip_system.update_occupation(get_checkerboard(&dip_system.latt));

        let interaction = get_dd_int_site(0, 0, &dip_system);
        assert_relative_eq!(interaction, 4.);  
    }

    #[test]
    fn m_matrix_det_test() {
        let mut dip_system = DipolarSystem::new(0., 0., 20., 1, 2);

        dip_system.update_occupation(get_checkerboard(&dip_system.latt));

        generate_dd_int_mat(&mut dip_system);
        
        let t = 1.;
        let mu = 1.;
        let m_mat = generate_mat_m(mu, t, &dip_system);
        
        // value taken from the python version of the code
        assert_relative_eq!(m_mat.determinant(), -0.4736842105263155);
    }

    #[test]
    fn get_tunneling_test() {
        let mut dip_system = DipolarSystem::new(0., 0., 20., 1, 2);
        dip_system.update_occupation(get_checkerboard(&dip_system.latt));
        generate_dd_int_mat(&mut dip_system);

        let mu = 1.;
        let tunneling = get_tunneling(mu, &dip_system, 1., 1e-3, 1e-2);

        // value taken from the python version of the code
        assert_relative_eq!(tunneling, 0.8200000000000006);

    }

    #[test]
    fn get_mu_inequality_test() {
        let mut dip_system = DipolarSystem::new(0., 0., 20., 1, 2);
        dip_system.update_occupation(get_checkerboard(&dip_system.latt));
        generate_dd_int_mat(&mut dip_system);

        let (lower, upper) = get_mu_inequality(&dip_system);
        assert_relative_eq!(lower, 0.);
        assert_relative_eq!(upper, 4.);
    }
}
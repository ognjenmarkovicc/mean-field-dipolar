#[macro_use]
extern crate approx; // For the macro relative_eq!
extern crate nalgebra as na;
use na::{Vector3};

/// Get the dipole-dipole interaction
/// 
/// dist_v is the distance vector
/// dip_v is the dipole vector, assumed to be a unit vector
pub fn get_dd_int(dist_v: Vector3<f64>, dip_v: Vector3<f64>)
-> f64 {
    let dist = dist_v.norm();

    let dd_int = (1.-3.*(dist_v.dot(&dip_v)/dist).powi(2))
        /(dist.powi(3));

    dd_int
}
/// Struct that holds periodic lattice information
pub struct PeriodicLattice 
{
    pub system_size: i32,
}

impl PeriodicLattice {
    /// Get a periodic lattice index
    /// from a bare lattice index
    /// 
    /// # Examples
    /// ```
    /// use mean_field_dipolar::PeriodicLattice;
    /// let system = PeriodicLattice { system_size: 4 };
    /// let idx = 5;
    /// assert_eq!(1, system.get_idx_periodic(idx));
    /// ```
    pub fn get_idx_periodic(&self, idx: i32) -> i32 {
        idx%self.system_size
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn periodic_latt_test() {
        let system = PeriodicLattice { system_size: 4 };
        let idx = 5;
        assert_eq!(1, system.get_idx_periodic(idx));
    }
}
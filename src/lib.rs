#[macro_use]
extern crate approx; // For the macro relative_eq!
extern crate nalgebra as na;
use na::{Vector3, DMatrix};
use std::f64::consts::PI;

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
/// Struct that holds periodic 2d lattice information
#[derive(Debug)]
pub struct PeriodicLattice 
{
    pub system_size: i32, // the system is a square
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
        idx.rem_euclid(self.system_size)
    }
}

/// Struct holding info about the dipolar system parameters
pub struct DipolarSystem {
    pub theta: f64,
    pub phi: f64,
    pub u_onsite: f64, // onsite interaction
}

impl DipolarSystem {
    pub fn get_dipole_vec(&self) -> Vector3<f64> {
        Vector3::new(self.theta.sin()*self.phi.cos(),
                     self.theta.sin()*self.phi.sin(),
                     self.theta.cos())
    }
}

/// Struct representing the pair of lattice positions
#[derive(Debug)]
pub struct LattPos<'a> {
    x: i32, // j
    y: i32, // i
    latt: &'a PeriodicLattice,
}

impl <'a> LattPos<'a> {
    /// Create a new lattice position
    /// 
    /// Panics if the given lattice indices are not compatible with
    /// latt.system_size.
    pub fn new(x: i32, y: i32, latt: &PeriodicLattice) -> LattPos {
        if x >= latt.system_size || y >= latt.system_size || x < 0 || y < 0 {
            panic!("Given indices not compatible with the given system size.");
        }
        LattPos { x: x, y: y, latt: latt }
    }
}

/// Struct representing the spin index 
/// of a lattice site
pub struct SpinIdx<'a> {
    idx: i32,
    latt: &'a PeriodicLattice,
}

impl <'a> SpinIdx<'a> {
    /// Create a new spin index
    /// 
    /// Panics if the index is not compatible with the given lattice.
    pub fn new(idx: i32, latt: &PeriodicLattice) -> SpinIdx {
        if idx >= latt.system_size.pow(2) || idx < 0 {
            panic!("Index not compatible with the given system size");
        }
        SpinIdx {idx, latt}
    }
}

impl <'a> From<SpinIdx<'a>> for LattPos<'a> {
    /// Get a lattice position from a spin index
    fn from(sp: SpinIdx) -> LattPos {
        LattPos { x: sp.idx%sp.latt.system_size,
                  y: sp.idx/sp.latt.system_size,
                  latt: sp.latt }
    }
}

impl <'a> From<LattPos<'a>> for SpinIdx<'a> {
    /// Get a spin index from a position in the lattice
    fn from(pos: LattPos) -> SpinIdx {
        SpinIdx {idx: pos.y*pos.latt.system_size + pos.x,
                 latt: pos.latt}
    }
}

pub fn get_checkerboard(latt: &PeriodicLattice) -> DMatrix<u16> {
    let mut mat: DMatrix<u16> = DMatrix::zeros(latt.system_size.try_into().unwrap(),
                                    latt.system_size.try_into().unwrap());

    for i in 0..latt.system_size {
        for j in 0..latt.system_size {
            mat[(i.try_into().unwrap(), j.try_into().unwrap())] = ((i + j)%2).try_into().unwrap();
        }
    }

    mat
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
    fn periodic_idx_pos_test() {
        let system = PeriodicLattice { system_size: 4 };
        let idx = 5;
        assert_eq!(1, system.get_idx_periodic(idx));
    }

    #[test]
    fn periodic_idx_neg_test() {
        let system = PeriodicLattice { system_size: 4 };
        let idx = -1;
        assert_eq!(3, system.get_idx_periodic(idx));
    }

    #[test]
    #[should_panic]
    fn latt_pos_panic_test() {
        let system = PeriodicLattice { system_size: 4 };
        let pos = LattPos::new(6, 5, &system);
    }

    #[test]
    fn pos_to_spin() {
        let system = PeriodicLattice { system_size: 4 };
        let pos = LattPos::new(3, 1, &system);
        let sp = SpinIdx::from(pos);

        assert_eq!(sp.idx, 7);
    }

    #[test]
    fn spin_to_pos() {
        let system = PeriodicLattice { system_size: 4 };
        let sp = SpinIdx::new(7, &system);
        let pos = LattPos::from(sp);

        assert_eq!(pos.x, 3);
        assert_eq!(pos.y, 1);
    }

    #[test]
    fn dipole_vec_test() {
        let dip_system = DipolarSystem {theta: PI/2., phi: 0., u_onsite: 0.};
        let dip_vec = dip_system.get_dipole_vec();

        relative_eq!(dip_vec[0], 1.);
    }
}
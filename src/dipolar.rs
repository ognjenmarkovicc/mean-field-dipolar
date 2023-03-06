use na::{Vector3, DMatrix};

use crate::lattice::PeriodicLattice;

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

/// Struct holding info about the dipolar system parameters
pub struct DipolarSystem {
    pub theta: f64,
    pub phi: f64,
    pub u_onsite: f64, // onsite interaction
    pub int_range: usize,
}

impl DipolarSystem {
    pub fn get_dipole_vec(&self) -> Vector3<f64> {
        Vector3::new(self.theta.sin()*self.phi.cos(),
                     self.theta.sin()*self.phi.sin(),
                     self.theta.cos())
    }
}

/// Get the dipole dipole interaction
/// a particle would experience if added to site (y, x)
pub fn get_dd_int_site(x: isize, y: isize,
                       dip: &DipolarSystem,
                       occupation: &DMatrix<u8>,
                       latt: &PeriodicLattice) -> f64 {

    let mut interaction: f64 = 0.;
    let int_range_cast = isize::try_from(dip.int_range).unwrap();

    // go over neighbors
    for x_n in x-int_range_cast..x+int_range_cast + 1 {
        for y_n in y-int_range_cast..y+int_range_cast + 1 {
            // get periodic indices
            let x_n_p = latt.get_idx_periodic(x_n);
            let y_n_p = latt.get_idx_periodic(y_n);

            let dist_vec =  Vector3::new(x as f64-x_n as f64,
                                         y as f64 -y_n as f64,
                                         0.);

            let dist = dist_vec.norm();

            if !(x_n == x && y_n == y) && dist<=dip.int_range as f64 {
                interaction += occupation[(y_n_p, x_n_p)] as f64 * get_dd_int(dist_vec, dip.get_dipole_vec());
            }
        }
    }

    interaction
}
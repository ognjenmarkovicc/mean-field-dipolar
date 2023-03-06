use na::{Vector3, DMatrix};
use nalgebra::{RowDVector, DVector};

use crate::lattice::{PeriodicLattice, SpinIdx, LattPos};

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
/// Assume that occupation is a latt.system_size x latt.system_size
/// matrix.
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

            if !(x_n_p == x as usize && y_n_p == y as usize) && dist<=dip.int_range as f64 {
                interaction += occupation[(y_n_p, x_n_p)] as f64 * get_dd_int(dist_vec, dip.get_dipole_vec());
            }
        }
    }

    interaction
}

// Generate the d-d interaction matrix
/// Assume that occupation is a latt.system_size x latt.system_size
/// matrix.
pub fn generate_dd_int_mat(dip: &DipolarSystem,
                           occupation: &DMatrix<u8>,
                           latt: &PeriodicLattice) -> DMatrix<f64> {

    let mut dd_mat: DMatrix<f64> 
        = DMatrix::zeros(latt.system_size, latt.system_size);
    for x in 0..latt.system_size {
        for y in 0..latt.system_size {
            // x, y guarnateed to be isize as
            // latt.system_size is < isize::MAX
            dd_mat[(y, x)] = get_dd_int_site(x as isize, y as isize, dip, occupation, latt);
        }
    }

    dd_mat
}

pub fn get_particle_e(x: usize, y: usize, mu: f64, dip: &DipolarSystem,
                      occupation: &DMatrix<u8>,
                      dd_mat: &DMatrix<f64>) -> f64 {
    -mu + dip.u_onsite*(occupation[(y, x)] as f64) + dd_mat[(y, x)]
}

pub fn get_hole_e(x: usize, y: usize, mu: f64, dip: &DipolarSystem,
    occupation: &DMatrix<u8>,
    dd_mat: &DMatrix<f64>) -> f64 { 
    mu - dip.u_onsite*(occupation[(y, x)] as f64 - 1.) - dd_mat[(y, x)]
}

pub fn get_m_row(spin_idx: &SpinIdx, mu: f64, t: f64,
                 dip: &DipolarSystem,
                 occupation: &DMatrix<u8>,
                 latt: &PeriodicLattice,
                 dd_mat: &DMatrix<f64>) -> RowDVector<f64> {

    let mut mat_row: RowDVector<f64>
        = RowDVector::zeros(latt.system_size*latt.system_size);
    mat_row[spin_idx.idx] = 1.;

    let latt_pos = LattPos::from(spin_idx);
    let n = occupation[(latt_pos.y, latt_pos.x)] as f64;

    let x = latt_pos.x as isize;
    let y = latt_pos.y as isize;

    let particle_e = get_particle_e(latt_pos.x, latt_pos.y,
                                    mu, dip, occupation, dd_mat);
    let hole_e = get_hole_e(latt_pos.x, latt_pos.y,
                            mu, dip, occupation, dd_mat);

    let row_val = if particle_e==0. || hole_e==0. {
        f64::INFINITY
    } else {
        -t*((n + 1.)/particle_e + n/hole_e) 
    };

    for (x_n, y_n) in [(x, y-1), (x, y+1), (x-1, y), (x+1, y)] {
        
        // get periodic indices
        let x_n_p = latt.get_idx_periodic(x_n);
        let y_n_p = latt.get_idx_periodic(y_n);

        // get the spin index of the neighbor
        let spin_idx_n = SpinIdx::from(LattPos::new(x_n_p, y_n_p, &latt)); 
        mat_row[spin_idx_n.idx] = row_val
    }

    mat_row
}

pub fn generate_mat_m(mu: f64, t: f64,
                      dip: &DipolarSystem,
                      occupation: &DMatrix<u8>,
                      latt: &PeriodicLattice,
                      dd_mat: &DMatrix<f64>) -> DMatrix<f64> {

    let mut m_mat = DMatrix::zeros(latt.system_size.pow(2),
                                   latt.system_size.pow(2));

    for spin_idx in 0..latt.system_size.pow(2) {

        let mat_row = get_m_row(&SpinIdx::new(spin_idx, latt), mu, t,
                                dip, occupation, latt, dd_mat);

        m_mat.set_row(spin_idx, &mat_row);
    }

    m_mat
}


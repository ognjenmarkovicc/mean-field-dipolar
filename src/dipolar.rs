use std::{path::{Path, PathBuf}};

use na::{Vector3, DMatrix};
use nalgebra::{RowDVector, DVector};

use crate::lattice::{PeriodicLattice, SpinIdx, LattPos, get_checkerboard};
use crate::util;

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
#[non_exhaustive]
pub struct DipolarSystem {
    pub theta: f64,
    pub phi: f64,
    pub u_onsite: f64, // onsite interaction
    pub int_range: usize,
    pub latt: PeriodicLattice,
    pub occupation: DMatrix<u8>,
    pub dd_mat: DMatrix<f64>,
}

impl DipolarSystem {
    pub fn get_dipole_vec(&self) -> Vector3<f64> {
        Vector3::new(self.theta.sin()*self.phi.cos(),
                     self.theta.sin()*self.phi.sin(),
                     self.theta.cos())
    }

    pub fn new(theta: f64, phi: f64, u_onsite: f64, 
               int_range: usize, system_size: usize
               ) -> Self {

        let latt = PeriodicLattice::new(system_size);

        let occupation = DMatrix::zeros(system_size, system_size);
        let dd_mat = DMatrix::zeros(system_size, system_size);
        DipolarSystem { theta, phi, u_onsite, int_range, latt, occupation, dd_mat}
    }

    pub fn update_occupation(&mut self, occupation: DMatrix<u8>) {
        assert!(occupation.nrows() == self.latt.system_size, "occupation nrows != system_size");
        
        assert!(occupation.ncols() == self.latt.system_size, "occupation ncols != system_size");

        self.occupation = occupation
    }
}

/// Get the dipole dipole interaction
/// a particle would experience if added to site (y, x)
/// Assume that occupation is a latt.system_size x latt.system_size
/// matrix.
pub fn get_dd_int_site(x: isize, y: isize,
                       dip: &DipolarSystem) -> f64 {

    let mut interaction: f64 = 0.;
    let int_range_cast = isize::try_from(dip.int_range).unwrap();

    // go over neighbors
    for x_n in x-int_range_cast..x+int_range_cast + 1 {
        for y_n in y-int_range_cast..y+int_range_cast + 1 {
            // get periodic indices
            let x_n_p = dip.latt.get_idx_periodic(x_n);
            let y_n_p = dip.latt.get_idx_periodic(y_n);

            let dist_vec =  Vector3::new(x as f64-x_n as f64,
                                         y as f64 -y_n as f64,
                                         0.);

            let dist = dist_vec.norm();

            if !(x_n_p == x as usize && y_n_p == y as usize) && dist<=dip.int_range as f64 {
                interaction += dip.occupation[(y_n_p, x_n_p)] as f64 * get_dd_int(dist_vec, dip.get_dipole_vec());
            }
        }
    }

    interaction
}

// Generate the d-d interaction matrix
/// Assume that occupation is a latt.system_size x latt.system_size
/// matrix.
pub fn generate_dd_int_mat(dip: &mut DipolarSystem) {

    let l = dip.latt.system_size;
    let mut dd_mat: DMatrix<f64> 
        = DMatrix::zeros(l, l);
    for x in 0..l {
        for y in 0..l {
            // x, y guarnateed to be isize as
            // l is < isize::MAX
            dd_mat[(y, x)] = get_dd_int_site(x as isize, y as isize, dip);
        }
    }

    dip.dd_mat = dd_mat;
}

pub fn get_particle_e(x: usize, y: usize, mu: f64, dip: &DipolarSystem) -> f64 {
    -mu + dip.u_onsite*(dip.occupation[(y, x)] as f64) + dip.dd_mat[(y, x)]
}

pub fn get_hole_e(x: usize, y: usize, mu: f64, dip: &DipolarSystem) -> f64 { 
    mu - dip.u_onsite*(dip.occupation[(y, x)] as f64 - 1.) - dip.dd_mat[(y, x)]
}

/// Generate Matrix M row from Trefzger et al., J. Phys. B At. Mol. Opt. Phys. 44 (2011) 193001, Eq. 3.19
/// 
/// # Parameters:
/// * mu - chemical potential
/// * t - tunneling
pub fn get_m_row(spin_idx: &SpinIdx, mu: f64, t: f64,
                 dip: &DipolarSystem) -> RowDVector<f64> {

    let latt = &dip.latt;
    let mut mat_row: RowDVector<f64>
        = RowDVector::zeros(latt.system_size*latt.system_size);
    mat_row[spin_idx.idx] = 1.;

    let latt_pos = LattPos::from(spin_idx);
    let n = dip.occupation[(latt_pos.y, latt_pos.x)] as f64;

    let x = latt_pos.x as isize;
    let y = latt_pos.y as isize;

    let particle_e = get_particle_e(latt_pos.x, latt_pos.y, mu, dip);
    let hole_e = get_hole_e(latt_pos.x, latt_pos.y, mu, dip);

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

/// Generate Matrix M from Trefzger et al., J. Phys. B At. Mol. Opt. Phys. 44 (2011) 193001, Eq. 3.19
/// 
/// # Parameters:
/// * mu - chemical potential
/// * t - tunneling
pub fn generate_mat_m(mu: f64, t: f64,
                      dip: &DipolarSystem) -> DMatrix<f64> {

    let latt = &dip.latt;
    let mut m_mat = DMatrix::zeros(latt.system_size.pow(2),
                                   latt.system_size.pow(2));

    for spin_idx in 0..latt.system_size.pow(2) {

        let mat_row = get_m_row(&SpinIdx::new(spin_idx, latt), mu, t, dip);

        m_mat.set_row(spin_idx, &mat_row);
    }

    m_mat
}

/// Find smallest tunneling where the det(M)=0,
/// where matrix M is defined in: 
/// Trefzger et al., J. Phys. B At. Mol. Opt. Phys. 44 (2011) 193001
/// Eq. 3.19
pub fn get_tunneling(mu: f64, dip: &DipolarSystem, max_tunneling: f64,
    tunneling_step:f64, det_threshold:f64) -> f64 {

    let mut tunneling = 1e-2;

    while tunneling < max_tunneling {
        let m_mat = generate_mat_m(mu, tunneling, dip);
        let det_val = m_mat.determinant();

        if det_val.abs() < det_threshold {
            return tunneling
        }
        tunneling += tunneling_step;
    }

    // return 0 if determinant smaller than det_threshold not found
    0.
}

/// Get the range of mu values in which the 
/// occupation is stable with the parameters
/// given in dip
pub fn get_mu_inequality(dip: &DipolarSystem) -> (f64, f64) {

    let n_float = dip.occupation.clone().cast::<f64>();
    let lower = dip.u_onsite*(n_float.add_scalar(-1.)) + &dip.dd_mat;
    let upper = dip.u_onsite*n_float + &dip.dd_mat;

    (lower.max(), upper.min())
}

pub fn simulation_sweep<P: AsRef<Path>>(save_path: P, int_ranges: (usize, usize), 
                        system_sizes: (usize, usize), theta: f64, phi: f64,
                        u_onsite: f64) {

    for int_range in (int_ranges.0..int_ranges.1).step_by(1) {
        for system_size in (system_sizes.0..system_sizes.1).step_by(2) {
            println!("Running int range {}, system size {}", int_range, system_size);

            let mut dip_system = DipolarSystem::new(theta, phi, u_onsite, int_range, system_size);
            dip_system.update_occupation(get_checkerboard(&dip_system.latt));
            generate_dd_int_mat(&mut dip_system);
            let (lower, upper) = get_mu_inequality(&dip_system);

            println!("Lower mu {:.2} upper mu {:.2}", lower, upper);

                if lower < upper {
                    let no_points = 100;
                    let mu_vals = util::linspace(lower, upper, no_points, true);

                    let tunneling = DVector::from_iterator(no_points, 
                                                        mu_vals.iter()
                                                                .map(|mu| get_tunneling(*mu, &dip_system,
                                                                                        4.,  1e-3, 1e-2)));
    
                    util::save_vector_json(save_path.as_ref()
                                                    .join(format!("tunneling_{system_size}_range_{int_range}.json"))
                                            , tunneling);
                    util::save_vector_json(save_path.as_ref()
                                                    .join(format!("mu_{system_size}_range_{int_range}.json"))
                                           , mu_vals);
            }
        }
    }
}
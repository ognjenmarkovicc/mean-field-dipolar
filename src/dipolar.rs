use na::Vector3;

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
    pub interaction_range: u16,
}

impl DipolarSystem {
    pub fn get_dipole_vec(&self) -> Vector3<f64> {
        Vector3::new(self.theta.sin()*self.phi.cos(),
                     self.theta.sin()*self.phi.sin(),
                     self.theta.cos())
    }
}
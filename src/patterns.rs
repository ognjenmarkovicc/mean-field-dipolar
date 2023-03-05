use na::{Vector3, DMatrix};
use super::lattice::PeriodicLattice;

pub fn get_checkerboard(latt: &PeriodicLattice) -> DMatrix<u8> {
    let mut mat: DMatrix<u8> = DMatrix::zeros(latt.system_size,
                                              latt.system_size);

    for i in 0..latt.system_size {
        for j in 0..latt.system_size {
            // assume usize is at least u32
            mat[(i, j)] = ((i + j)%2).try_into().unwrap();
        }
    }

    mat
}
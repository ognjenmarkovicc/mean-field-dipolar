use na::{DMatrix};
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

pub fn get_filled(latt: &PeriodicLattice) -> DMatrix<u8> {
    DMatrix::repeat(latt.system_size, latt.system_size, 1)
}

pub fn get_horizontal_stripe(latt: &PeriodicLattice) -> DMatrix<u8> {
    DMatrix::from_fn(latt.system_size, latt.system_size, |row, _| (row%2) as u8)
}
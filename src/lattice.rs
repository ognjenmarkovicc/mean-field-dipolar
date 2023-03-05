

/// Struct that holds periodic 2d lattice information
#[derive(Debug)]
#[non_exhaustive]
pub struct PeriodicLattice 
{
    pub system_size: usize, // the system is a square
}

impl PeriodicLattice {
    /// Get a periodic lattice index
    /// from a bare lattice index
    /// 
    /// # Examples
    /// ```
    /// use mean_field_dipolar::lattice::PeriodicLattice;
    /// let system = PeriodicLattice::new(4);
    /// let idx = 5;
    /// assert_eq!(1, system.get_idx_periodic(idx));
    /// ```
    pub fn get_idx_periodic(&self, idx: isize) -> usize {
        // can get a negative argument
        // convert system size into isize first and then
        // convert the result to usize
        (idx.rem_euclid(self.system_size.try_into().unwrap()))
            .try_into().unwrap()
    }

    pub fn new(system_size: usize) -> Self {
        if system_size > isize::MAX.try_into().unwrap() {
            panic!("Given system_size needs to be less than {}", isize::MAX);
        }

        PeriodicLattice { system_size }
    }
}

/// Struct representing the pair of lattice positions
#[derive(Debug)]
#[non_exhaustive]
pub struct LattPos<'a> {
    pub x: usize, // j
    pub y: usize, // i
    pub latt: &'a PeriodicLattice,
}

impl <'a> LattPos<'a> {
    /// Create a new lattice position
    /// 
    /// Panics if the given lattice indices are not compatible with
    /// latt.system_size.
    pub fn new(x: usize, y: usize, latt: &PeriodicLattice) -> LattPos {
        if x >= latt.system_size || y >= latt.system_size || x < 0 || y < 0 {
            panic!("Given indices not compatible with the given system size.");
        }
        LattPos { x: x, y: y, latt: latt }
    }
}

/// Struct representing the spin index 
/// of a lattice site
#[non_exhaustive]
pub struct SpinIdx<'a> {
    pub idx: usize,
    pub latt: &'a PeriodicLattice,
}

impl <'a> SpinIdx<'a> {
    /// Create a new spin index
    /// 
    /// Panics if the index is not compatible with the given lattice.
    pub fn new(idx: usize, latt: &PeriodicLattice) -> SpinIdx {
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
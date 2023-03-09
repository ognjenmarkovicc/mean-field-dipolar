use na::{DVector};

/// Basic linspace function
/// 
/// Doesn't perform any checks, use with
/// care. Will silently convert num and idx to f64 
/// even if num cannot be represented as f64
pub fn linspace(start: f64, stop: f64, num: usize, endpoint: bool) -> DVector<f64> {

    let denom = if endpoint {
        num - 1
    } else {
        num
    };

    let delta = (stop - start)/(denom as f64);
    DVector::from_iterator(num, (0..num).map(|idx| (idx as f64)*delta + start))
}
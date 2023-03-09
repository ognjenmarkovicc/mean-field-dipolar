use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mean_field_dipolar::dipolar::{DipolarSystem, generate_dd_int_mat, generate_mat_m};
use mean_field_dipolar::patterns::get_checkerboard;

pub fn mmat_det_benchmark_fn(system_size: usize) {
    let mut dip_system = DipolarSystem::new(0., 0., 20., 1, system_size);
    dip_system.update_occupation(get_checkerboard(&dip_system.latt));
    generate_dd_int_mat(&mut dip_system);

    let m_mat = generate_mat_m(1., 1., &dip_system);
    let det_val = m_mat.determinant();
}

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("mmat det", |b| b.iter(|| mmat_det_benchmark_fn(black_box(4))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

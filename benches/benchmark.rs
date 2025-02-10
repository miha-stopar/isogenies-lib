use criterion::{criterion_group, criterion_main, Criterion};
use isogenies_lib::schemes::klapoti_test::{params128, params32, params512, params64};

pub fn criterion_benchmark(c: &mut Criterion) {
    // c.bench_function("params 32", |b| b.iter(|| params32()));
    // c.bench_function("params 64", |b| b.iter(|| params64()));
    // c.bench_function("params 128", |b| b.iter(|| params128()));
    // c.bench_function("params 512", |b| b.iter(|| params512()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

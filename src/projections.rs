// Lie Algebra sum and products of the Eij matrices
type LASum<T> = Vec<T>;
type LAProd<T> = Vec<T>;

fn find_casimir(p: usize, k: usize) -> LASum<LAProd<(usize, usize)>> {
    let mut q = Vec::new();
    q.push(Vec::new());
    let mut cycles = Vec::new();
    while !q.is_empty() {
        let str = q.remove(0);
        for i in 0..k {
            let mut str_clone = str.clone();
            str_clone.push(i);
            if str_clone.len() == p {
                cycles.push(str_clone);
            } else {
                q.push(str_clone);
            }
        }
    }
    let mut summation = Vec::new();
    for cycle in cycles {
        let mut entry = Vec::new();
        for i in 1..cycle.len() {
            entry.push( (cycle[i-1], cycle[i]) );
        }
        entry.push( (cycle[cycle.len()-1], cycle[0]) );
        summation.push(entry);
    }
    summation
}

use nalgebra::DMatrix;
fn casimir_eigenval(casimir_num: u32, lambda: &Vec<u32>, k: usize) -> f64 {
    let a = DMatrix::<f64>::from_fn(k, k, |i, j| {
        if i==j {
            lambda[i] as f64 + k as f64 - i as f64
        } else if i<j {
            -1.
        } else {
            0.
        }
    });
    let a_pow_p = a.pow(casimir_num);
    let e = DMatrix::<f64>::from_element(k, k, 1.);

    (a_pow_p * e).trace()
}

fn all_partitions(n: usize, m: usize) -> Vec<Vec<u32>> {
    // P1
    let mut part = vec![1;m+1];
    part[0] = 0;
    // P2
    part[m] = n;
    let mut q = m - (if n==1 {1} else {0});
    // P3
    



}
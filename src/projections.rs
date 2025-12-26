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

fn find_distinguishing_casimir_index(new_partition: &Vec<u32>, proj_partition: &Vec<u32>) -> Option<usize> {
    let k = new_partition.len();
    for p in 1..=k {
        if casimir_eigenval(p as u32, new_partition, k) != casimir_eigenval(p as u32, proj_partition, k) {
            return Some(p);
        }
    }
    None
}

fn find_distinguishing_parts_and_indices(proj_part: &Vec<u32>) -> Vec<(Vec<u32>, usize)> {
    let k = proj_part.len();
    let n: usize = proj_part.iter().sum::<u32>().try_into().unwrap();

    partitions_n_le_m(n, k).iter().filter(|p| find_distinguishing_casimir_index(p, proj_part).is_some())
                                    .map(|p| (p.clone(), find_distinguishing_casimir_index(p, proj_part).unwrap()))
                                    .collect()
}

// algorithm from pg 392 of knuth art of computer programming, combinatorial algorithms part 1
fn partitions_n_m(n: usize, m: usize) -> Vec<Vec<u32>> {
    let mut part = vec![1;m+1];
    let mut parts: Vec<Vec<u32>> = vec![];
    part[0] = n - m + 1;
    part[m] = 0;

    let mut j: usize;
    let mut s: usize;
    let mut x: usize;

    loop {
        parts.push( part.iter().map(|x| (*x).try_into().unwrap() ).collect() );
        // H2
        if part[1] < part[0] - 1 {
            // H3
            part[0] -= 1;
            part[1] += 1;
            continue;
        }
        // H4
        j=2;
        s=part[0]+part[1]-1;
        while part[j] >= part[0] - 1 {
            s += part[j];
            j += 1;
        }
        // H5
        if j >= m { break; }
        x = part[j] + 1;
        part[j] = x;
        j -= 1;
        // H6
        while j > 0 {
            part[j] = x;
            s -= x;
            j -= 1; 
        }
        part[0] = s;
    }
    parts.iter_mut().for_each(|p| { p.pop(); });
    parts
}

fn partitions_n_le_m( n: usize, m: usize) -> Vec<Vec<u32>> {
    let mut parts: Vec<Vec<u32>> = vec![ vec![0;m]];
    parts[0][0] = n.try_into().unwrap();
    let parts_len_lm_iterator = (2..=m).map(|lm| partitions_n_m(n, lm).into_iter().map(|mut p| {p.resize(m,0); p} ) );
    for parts_len_lm in parts_len_lm_iterator {
        parts.extend(parts_len_lm.into_iter());
    }
    parts
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn find_casimir_2_2() {
        let casimir = find_casimir(2, 2);
        print!(" casimir: {casimir:?}");
        assert!(true);
    }

    #[test]
    fn find_all_6_3_partitions_with_dist_casimir() {
        for part in partitions_n_le_m(6, 3) {
            if find_distinguishing_casimir_index(&part, &vec![6,0,0]).is_some() {
                println!("{part:?}");
            }
        }
        assert!(true);
    }

    #[test]
    fn find_partitions_11_4() {
        let parts = partitions_n_m(11, 4);
        print!(" parts: {parts:?}");
        assert!(true);
    }

    #[test]
    fn find_partitions_11_le_4() {
        let parts = partitions_n_le_m(11, 4);
        print!(" parts: {parts:?}");
        assert!(true);
    }
}

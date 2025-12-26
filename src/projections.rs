use super::straight_line_program::SLP;
use super::straight_line_program::{add_slp, mult_slp, scale_slp};
use super::straight_line_program::stringify_slp;

use super::transformations::transform_x_i_product_program;

// Lie Algebra sum and products of the Eij matrices
type LASum<T> = Vec<T>;
type LAProd<T> = Vec<T>;

fn find_casimir(p: usize, k: usize) -> LASum<LAProd<(usize, usize)>> {
    let mut q = Vec::new();
    q.push(Vec::new());
    let mut cycles = Vec::new();
    while !q.is_empty() {
        let str = q.remove(0);
        for i in 0..k-1 {
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
use num_complex::Complex64;
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

pub fn apply_projection_to_slp(slp: SLP, projection: &Vec<u32>) -> Option<SLP> {
    let distinguishing_partitions = find_distinguishing_parts_and_indices(projection);
    distinguishing_partitions.iter().for_each(|(part, el_index)| {
        println!("Found part {part:?}, with casimir el index: {el_index}");
    });
    
    let k = projection.len();
    let mut slp_product: Vec<SLP> = vec![];
    for (part, el_index) in distinguishing_partitions {
        let mut slp_term1 = slp.clone();    
        let mut slp_term2 = slp.clone();

        // apply U_{mu}
        let eij_sum_of_products = find_casimir(el_index, k);
        println!("Casimir found: {eij_sum_of_products:?}, index: {el_index}");
        let mut slp_term1_entries: Vec<SLP> = vec![];
        for eij_prod in  eij_sum_of_products {
            let transformed: Option<SLP> = transform_x_i_product_program(&slp_term1, &eij_prod, k);
            match transformed {
                Some(transform_slp) => slp_term1_entries.push(transform_slp),
                None => {
                    println!("Failure for projection {part:?} SLP:\n{}", stringify_slp(&slp_term1));
                }
            }
        }
        if slp_term1_entries.len() == 0 {
            continue;
        }
        slp_term1 = slp_term1_entries[0].clone();
        for i in 1..slp_term1_entries.len() {
            add_slp(&mut slp_term1, &mut slp_term1_entries[i]);
        }

        // transform term2 by chi_{mu}
        scale_slp(&mut slp_term2, 
            Complex64::new(-1.0*casimir_eigenval(el_index as u32, &part, k) , 0.0)
        );

        add_slp(&mut slp_term1, &mut slp_term2);
        scale_slp(&mut slp_term1, Complex64::new(
            casimir_eigenval(el_index as u32, projection, k)
             - casimir_eigenval(el_index as u32, &part, k)
             , 0.0
        ));

        slp_product.push(slp_term1);
    }

    let mut slp_proj = slp_product[0].clone();
    for i in 1..slp_product.len() {
        mult_slp(&mut slp_proj, &mut slp_product[i]);
    }
    Some(slp_proj)
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

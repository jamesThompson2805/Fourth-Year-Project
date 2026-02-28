//! projections.rs contains functions to achieve the projection of a slp to a lambda-isotypic space
#![allow(dead_code)]
use std::fmt::Display;

use super::straight_line_program::SLP;
use super::straight_line_program::{add_slp, scale_slp};
use super::straight_line_program::stringify_slp;
use super::pbw_reduce::*;

use crate::transformations::{self, *};
use std::fmt::Debug;
use std::ops::{Add, Sub, Mul, Div};
use std::hash::Hash;

/// LASum wraps a vector denoting a collection of terms to be summed together
type LASum<T> = Vec<T>;
/// LAProd wraps a vector denoting a collection of terms to be multiplied together
type LAProd<T> = Vec<T>;

/// find casimir finds the pth casimir element for the casimir elements of k variables
/// p is 1 indexed
fn find_casimir(p: usize, k: usize) -> LASum<LAProd<(usize, usize)>> {
    let mut q = Vec::new();
    q.push(Vec::new());
    let mut cycles = Vec::new();
    while !q.is_empty() {
        let str = q.remove(0);
        for i in 0..=k-1 {
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

use itertools::Itertools;
use nalgebra::DMatrix;
/// casimir_eigenval finds an eigenvalue of a casimir element wrt partition lambda, given by p, the index of the casimir number
/// casimir_num is the index of the casimir number
/// lambda is the partition to find the eigenvalue with respect to
/// k is the number of variables of the metapolynomial, and the length of the partition and the 
fn casimir_eigenval(casimir_num: u32, lambda: &Vec<u32>, k: usize) -> i32 {
    let a = DMatrix::<i32>::from_fn(k, k, |i, j| {
        if i==j {
            lambda[i] as i32 + k as i32 - i as i32 - 1
        } else if i<j {
            -1
        } else {
            0
        }
    });
    let a_pow_p = a.pow(casimir_num);
    let e = DMatrix::<i32>::from_element(k, k, 1);

    (a_pow_p * e).trace()
}

/// find_distinguishing_casimir_index finds the index of a casimir element that has a different eigenvalue for the two partitions provided
/// new_partition is the first partition
/// proj_partitions is the other partition
/// returns an Option of a usize as it may be the case that there doesn't exist a distinguishing casimir element (occurring when we choose partitions that aren't highest weight)
fn find_distinguishing_casimir_index(new_partition: &Vec<u32>, proj_partition: &Vec<u32>) -> Option<usize> {
    let k = new_partition.len();
    for p in 1..=k {
        if casimir_eigenval(p as u32, new_partition, k) != casimir_eigenval(p as u32, proj_partition, k) {
            return Some(p);
        }
    }
    None
}

/// find_distinguishing_parts_and_indices attempts to finds all other partitions that can be distinguished by some casimir element
/// proj_part is a partition
/// returns a vector of every other partition of the same length and the index of the distinguishing casimir element
fn find_distinguishing_parts_and_indices(proj_part: &Vec<u32>) -> Vec<(Vec<u32>, usize)> {
    let k = proj_part.len();
    let n: usize = proj_part.iter().sum::<u32>().try_into().unwrap();

    partitions_n_le_m(n, k).iter().filter(|p| find_distinguishing_casimir_index(p, proj_part).is_some())
                                    .map(|p| (p.clone(), find_distinguishing_casimir_index(p, proj_part).unwrap()))
                                    .collect()
}

/// partitions_n_m finds all partitions of n into exactly m parts
///  algorithm from pg 392 of knuth art of computer programming, combinatorial algorithms part 1
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

/// partitions_n_le_m finds all partitions of n into m parts or fewer, each part is padded to a length of m
fn partitions_n_le_m( n: usize, m: usize) -> Vec<Vec<u32>> {
    let mut parts: Vec<Vec<u32>> = vec![ vec![0;m]];
    parts[0][0] = n.try_into().unwrap();
    let parts_len_lm_iterator = (2..=m).map(|lm| partitions_n_m(n, lm).into_iter().map(|mut p| {p.resize(m,0); p} ) );
    for parts_len_lm in parts_len_lm_iterator {
        parts.extend(parts_len_lm.into_iter());
    }
    parts
}

/// apply_casimir_to_slp_slow applies a casimir element to a slp and potentially returns the transformation (or an error message)
/// It works by applying each Eij transformation one after another, generating exponentially many lines!
/// 
/// slp is the slp
/// casimir_index is the index of the casimir element (starting at 1)
/// casimir_size should match the number of variables the metapolynomial slp would utilise
/// u32_to_c is a function converting a u32 to the custom type T of the SLP (should be sensible, u32_to_c(0) -> T::ZERO )
fn apply_casimir_to_slp_slow<T,F>(slp: &SLP<T>, casimir_index: usize, casimir_size: usize, u32_to_c: F) -> Result<SLP<T>, String>
where T: Clone + Display, F: Fn(u32) -> T,
{
    let casimir = find_casimir(casimir_index, casimir_size);

    let mut mu_c_slp_summands: LASum<SLP<T>> = vec![];
    for eij_prod in casimir.iter() {
        let transformed_res= apply_eijs_on_program(&slp, &eij_prod, |x| u32_to_c(x));
        match transformed_res {
            Ok(transform_slp) => mu_c_slp_summands.push(transform_slp),
            Err(s) => return Err( format!("Failure at casimir {casimir:?}, error{s}, SLP:\n{}", stringify_slp(&slp))),
        }
        //println!("Transformed slp: \n{}", stringify_slp(mu_c_slp_summands.last().unwrap()));
    }
    if mu_c_slp_summands.len() == 0 {
        return Err(String::from("casimir element empty, invalid index given"));
    }
    let first_el = mu_c_slp_summands[0].clone();

    Ok(
        mu_c_slp_summands.into_iter().skip(1).fold(first_el, |acc, slp| add_slp(acc,slp))
    )
}

/// apply_projection_to_slp_slow applies the lambda projection to a slp
/// the projection is defined as P_{\lambda} = \Pi_{\mu \in \Lambda, \mu \neq \lambda} {
///     \frac{ U_{\mu} - \Chi_{\mu}(U_{\mu}) }{ \Chi_{\lambda}(U_{\mu}) - \Chi_{\mu}(U_{\mu}) }
/// }
/// 
/// when evaluated on a metapolynomial, U_{\mu} is applied to the metapolynomial, a scale of the metapoly by \Chi_{\mu}(U_{\mu}) is subtracted
///  and the result is scaled by the denominator. All terms are then summed.
/// 
/// This method is again denoted slow for using the naive method and generating exponentiall many lines
/// 
/// slp is the slp to project
/// lambda is the corresponding integer partition to project the slp to
/// i64_to_c is a function converting integers into the coefficient object, it should be sensible wrt numerical operations and default, i64_to_c(0) = 0
/// unit should be the multiplicative identity of the vector space being used for T
/// returns a result of either the projected slp or an error message detailing the issue
pub fn apply_projection_to_slp_slow<T, F>(slp: SLP<T>, lambda: &Vec<u32>, i32_to_c: F, unit: T) -> Result<SLP<T>, String>
where
    T: Clone + Display + Default + Add<Output=T> +Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Debug + Eq,
    F: Fn(i32) -> T,
{
    // steps of finding the projection:
    // 1. Find all candidate partitions
    // 2. for each partition, calculate its term in the product
    //  2.1 calculate the projection under U_{\mu}
    //  2.2 calculate the scale by -1 * \Chi_{\mu}(U_{\mu}) 
    //  2.3 scale result by \Chi_{\lambda}(U_{\mu}) - \Chi_{\mu}(U_{\mu})
    // 3. multiply all these terms together and return the resulting slp

    let distinguishing_partitions = find_distinguishing_parts_and_indices(lambda);
    distinguishing_partitions.iter().for_each(|(partition, el_index)| {
        println!("Found partition {partition:?}, with casimir el index: {el_index}");
    });

    apply_candidate_partitions_to_slp_slow(slp, lambda, &distinguishing_partitions, i32_to_c, unit)
}

use std::collections::HashMap;
/// eval_proj_pairs reduces the product of projections in eq17 to a list of projections (denoted by casimir elements) to add together,
///  this is necessary as we do not have distributivity and hence all need to be expanded out and computed as a sum of applications and scaling.
/// 
/// v is a vector of pairs (usize, T) representing a term in the product of the form (C_n, \Chi_{n}(C_n)) representing (C_n - \Chi_{n}(C_n))
/// unit is the representative of 1 in T
/// returns the list of projections to apply and sum together
fn eval_proj_pairs<T>(v: &[(usize, T)], unit: T) -> HashMap<LAProd<usize>, T>
where T: Mul<Output=T> + Add<Output=T> + Clone + Debug,
{
    if v.len() == 0 { HashMap::new() }
    else if v.len() == 1 { HashMap::from([ (vec![v[0].0], unit.clone()),  (vec![], v[0].1.clone()) ]) }
    else {
        let before = eval_proj_pairs(&v[1..], unit);
        let mut res: HashMap<Vec<usize>, T> = HashMap::new();
        for (vec,val) in before {
            if let Some(coeff) = res.get_mut(&vec) {
                *coeff = coeff.clone() + v[0].1.clone() * val.clone();
            } else {
                res.insert(vec.clone(), v[0].1.clone()*val.clone());
            }
            
            let vec_new = [vec, vec![v[0].0]].concat();
            if let Some(coeff) = res.get_mut(&vec_new) {
                *coeff = coeff.clone() + val;
            } else {
                res.insert(vec_new, val);
            }
        }
        res
    }
}

/// casimir_prod_to_sorted_basis_matrices takes a list of references to casimir elements and returns a polynomial of sorted basis transforms
fn casimir_prod_to_sorted_basis_matrices(prod: &[usize], k: usize) -> HashMap<LAProd<(usize,usize)>, i32>
{
    let counts = prod.iter().map(|u| find_casimir(*u, k)).multi_cartesian_product().map(|v| v.concat()).counts(); // have distributed out product and collected counts
    counts.into_iter().map(|(k,v)|{
        let mut sorted = pbw_reduce(&k[..]);
        sorted.iter_mut().for_each(|(_,vs)| *vs *= v as i32); 
        sorted
    }).flatten().into_grouping_map().sum()
}

/// casimir_eq_to_sorted_basis_matrices takes a sum of scaled products of casimir elements and computes their polynomial of sorted basis elements
fn casimir_eq_to_sorted_basis_matrices(sum: &Vec<(Vec<usize>,i32)>, k: usize) -> HashMap<LAProd<(usize, usize)>, i32> {
    sum.into_iter().map(|(v,c)| 
        casimir_prod_to_sorted_basis_matrices(v, k).into_iter().map(move |p| (p.0,p.1 * c))
    ).flatten().into_grouping_map().sum()
}

/// eval_proj_pairs_to_sorted_basis takes a list representing the casimir element and the term to add (from eq17. in paper)
///  and returns the polynomial of basis transformations to make
fn eval_proj_pairs_to_sorted_basis(v: &[(usize, i32)], k: usize) -> HashMap<LAProd<(usize,usize)>, i32> {
    let casimir_poly = eval_proj_pairs::<i32>(v, 1);
    println!("casimir poly {casimir_poly:?}");
    casimir_eq_to_sorted_basis_matrices(&casimir_poly.into_iter().collect(), k)
}


/// apply_candidate_partitions_to_slp_slow applies all candidate partitions to a slp in the exponentially slow way
/// the projection is defined as P_{\lambda} = \Pi_{\mu \in \Lambda, \mu \neq \lambda} {
///     \frac{ U_{\mu} - \Chi_{\mu}(U_{\mu}) }{ \Chi_{\lambda}(U_{\mu}) - \Chi_{\mu}(U_{\mu}) }
/// }
/// 
/// when evaluated on a metapolynomial, U_{\mu} is applied to the metapolynomial, a scale of the metapoly by \Chi_{\mu}(U_{\mu}) is subtracted
///  and the result is scaled by the denominator. All terms are then summed.
/// 
/// slp is the slp to project
/// lambda is the corresponding integer partition to project the slp to
/// partitions is a subset of all partitions of same length and integer as lambda alongside an index to a casimir element that distinguishes the two
/// i64_to_c is a function converting integers into the coefficient object, it should be sensible wrt numerical operations and default, i64_to_c(0) = 0
/// unit should be the multiplicative identity of the vector space being used for T
/// returns a result of either the projected slp or an error message detailing the issue
fn apply_candidate_partitions_to_slp_slow<T,F>(slp: SLP<T>, lambda: &Vec<u32>, partitions: &Vec<(Vec<u32>, usize)>, i32_to_c: F, unit: T) -> Result<SLP<T>, String>
where 
    T: Clone + Display + Default + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Add<Output=T> + Debug + PartialEq,
    F: Fn(i32) -> T,
{
    // steps of finding the projection:
    // 1. Find all candidate partitions (provided by partitions)
    // 2. Expand out the product of the partition terms
    // 3. Calculate the series of projections for each term in expanded product and scale appropriately
    // 4. Sum the products and scale by the denominator product
    if partitions.len() == 0 { return Err(String::from("No partitions found")); }
    let k = lambda.len();
    if partitions.iter().any(|(v,_)| v.len() != k) { return Err(String::from("Partitions not all same length")); }

    let terms: Vec<(usize, T)> = partitions.iter().map(|(v,i)| ( *i, (T::default() - unit.clone()) * i32_to_c(casimir_eigenval((*i) as u32, v, k)) )).collect();
    let expanded = eval_proj_pairs(&terms.as_slice(), unit.clone());
    let _ = dbg!(&expanded);

    let mut summands = Vec::new();
    for (proj_s, c) in expanded {
        // println!("Proj: {proj_s:?}");
        let mut slp_term = proj_s.into_iter().try_fold(slp.clone(),
         |acc, proj| apply_casimir_to_slp_slow(&acc, proj, k, |x| i32_to_c(x as i32))
        )?;
        scale_slp(&mut slp_term, c);
        summands.push( slp_term );
    }
    if summands.len() == 0 { return Err("no partitions provided a projection".into()); }
    let first = summands.first().unwrap().clone();
    let mut slp_res = summands.into_iter().skip(1).fold(first, |acc, term| add_slp(acc, term));

    let denom_scale = partitions.iter().map(
        |(v,i)| i32_to_c(casimir_eigenval((*i) as u32, lambda, k)) - i32_to_c(casimir_eigenval((*i) as u32, v, k))
    ).reduce(|acc,e| acc * e).ok_or("partitions is empty")?;

    scale_slp( &mut slp_res,  unit / denom_scale );

    Ok(slp_res)
}

pub fn apply_lambda_projection_to_slp<T, F, const K: usize>(
    slp: &SLP<T>,
    lambda: &Vec<u32>,
    i64_to_c: F,
) -> Result<SLP<T>, String>
where
    T: Add<Output=T> + Mul<Output=T> + Clone + Eq + Hash + Display + Debug,
    F: Fn(i64) -> T + Copy,
{
    
    let k = (K as f64).sqrt() as usize;
    if lambda.len() != k {
        return Err("The partitio given should have length equal to sqrt(K)".into());
    }

    let distinguishing_partitions: Vec<_> = find_distinguishing_parts_and_indices(&lambda).into_iter().map(|(p,c)| (p, c as u32)).collect();
    let numerators: Vec<(usize, i32)> = distinguishing_partitions.iter().map(|(p, cas_num)| (*cas_num as usize, casimir_eigenval(*cas_num, p, k))).collect();
    let denominator: i64 = distinguishing_partitions.iter().map(|(p, cas_num)| casimir_eigenval(*cas_num, &lambda, k) as i64 - casimir_eigenval(*cas_num, p, k) as i64 ).product();

    let numerator_basis_transforms = eval_proj_pairs_to_sorted_basis(&numerators[..], k);
    let mut slp_res = transformations::apply_eij_poly_on_program::<T,F,K>(slp, &numerator_basis_transforms, i64_to_c)?;
    
    let scalar = i64_to_c( 1 / denominator );
    scale_slp(&mut slp_res, scalar);

    Ok(slp_res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Rational64;
    use std::fs::File;
    use std::io::Read;
    use crate::evaluation::{are_rational_slp_similar, stepwise_slp_to_poly};
    use crate::straight_line_program::scale_slp;
    use crate::parsing::slp_parser_rational;

    #[test]
    fn find_casimir_2_2() {
        let casimir = find_casimir(2, 3);
        print!(" casimir: {casimir:?}");
        assert!(true);
    }

    #[test]
    fn find_casimir_eigenval_3_3() {
        assert_eq!( 48, casimir_eigenval(2, &vec![6,0,0], 3));
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

    #[test]
    fn test_projection_example_4_14_1_slow() {
        let slp_str = "=C<1,1>
=C<0,2>
=C<2,0>
=L0*L0
=L1*L2
=2*L4
=L3+L5";
        let slp = slp_parser_rational(slp_str).unwrap();
        let u32_to_c = |x:u32| Rational64::new(x.into(),1);
        let slp_proj = apply_casimir_to_slp_slow(&slp, 1, 2, u32_to_c).unwrap();
        let slp_proj2 = apply_casimir_to_slp_slow(&slp, 2, 2, u32_to_c).unwrap();

        let mut slp_res = slp.clone();
        let mut slp_res2 = slp.clone();
        scale_slp(&mut slp_res, Rational64::from_integer(4));
        scale_slp(&mut slp_res2, Rational64::from_integer(20));

        println!("{}", stepwise_slp_to_poly( &slp_proj2, Rational64::ONE ));

        assert!(are_rational_slp_similar(&slp_proj, &slp_res).1);
        assert!(are_rational_slp_similar(&slp_proj2, &slp_res2).1);
    }

    #[test]
    fn test_projection_to_example_5_4_slow() {
        let slp_str = "=C<0,0,2>
=C<0,2,0>
=C<2,0,0>
=L0*L1
=L3*L2";
        let slp = slp_parser_rational(slp_str).unwrap();
        let distinguishing_parts = vec![
            (vec![4,2,0], 2), (vec![2,2,2], 2)
        ];
        let lambda = vec![6,0,0];

        let mut res_str: String = String::new();
        File::open("example_5_4_proj.txt").unwrap().read_to_string(&mut res_str).unwrap();
        let res_slp = slp_parser_rational(&res_str).unwrap();
        // scale_slp(&mut res_slp, Rational64::new(1,60) );

        let i32_to_c = |i| { Rational64::new(i as i64,1) };
        let unit = Rational64::ONE;
        let mut proj_slp = apply_candidate_partitions_to_slp_slow(slp, &lambda, &distinguishing_parts, i32_to_c, unit).unwrap();
        scale_slp(&mut proj_slp, Rational64::new(60,1) );
        // println!("{}", stepwise_slp_to_poly(&proj_slp, Rational64::ONE));
        println!("{}", stepwise_slp_to_poly(&res_slp, Rational64::ONE));
        println!("{:?}", stepwise_slp_to_poly(&proj_slp, Rational64::ONE).split("\n").last());
        assert!(are_rational_slp_similar(&res_slp, &proj_slp).1);
    }

    // This test is too slow to run!
    #[test]
    fn test_projection_example_5_4_is_same_as_distinguishing_partitions_slow() {
        let slp_str = "=C<0,0,2>
=C<0,2,0>
=C<2,0,0>
=L0*L1
=L3*L2";
        let slp = slp_parser_rational(slp_str).unwrap();
        let lambda = vec![6,0,0];
        let distinguishing_parts = find_distinguishing_parts_and_indices(&lambda);

        let i32_to_c = |i| { Rational64::new(i as i64,1) };
        let unit = Rational64::ONE;
        let mut proj_slp = apply_candidate_partitions_to_slp_slow(slp, &lambda, &distinguishing_parts, i32_to_c, unit).unwrap();
        scale_slp(&mut proj_slp, Rational64::new(60,1) );

        let mut res_str: String = String::new();
        File::open("example_5_4_proj.txt").unwrap().read_to_string(&mut res_str).unwrap();
        let res_slp = slp_parser_rational(&res_str).unwrap();

        println!("{:?}", stepwise_slp_to_poly(&proj_slp, Rational64::ONE).split("\n").last());
        assert!(are_rational_slp_similar(&res_slp, &proj_slp).1);
    }

    #[test]
    fn test_pbw_reduce_trivials() {
        let starting = &[(2,2),(1,1),(0,0)]; // needs to be sorted but commutator always 0
        let starting2 = &[(0,1),(1,2),(2,2)]; // already sorted
        assert_eq!(pbw_reduce(starting), HashMap::from([(vec![(0,0),(1,1),(2,2)],1)]) );
        assert_eq!(pbw_reduce(starting2), HashMap::from([(vec![(0,1),(1,2),(2,2)],1)]) );
    }

    #[test]
    fn test_pbw_reduce_nontrivials() {
        let starting = &[(2,1),(1,0),(1,2)]; 
        // [21,10,12] -> [10,21,12]+[20,12] -> [10,12,21]+[10,22]-[10,11]+[20,12] -> [10,12,21]+[10,22]-[10,11]+[12,20]-[10]
        println!("{:?}", pbw_reduce(starting));
    }

    #[test]
    fn test_casimir_prod_sorting() {
        let casimirs = vec![(vec![1,2,3],1), (vec![2,1],1)];
        println!("{:?}", casimir_eq_to_sorted_basis_matrices(&casimirs, 3).len());
    }

    #[test]
    fn test_apply_casimir_on_program() {
        let slp_str = "=C<2,0,0>
=C<0,2,0>
=C<0,0,2>
=L0*L1
=L3*L2";
        let slp = slp_parser_rational(slp_str).expect("Should parse SLP");
        let i64_to_c = |i| Rational64::new(i, 1);
        let casimir_sorted = eval_proj_pairs_to_sorted_basis(&vec![(2,-12),(2,-28)], 3);

        println!("Casimir el 2 is {casimir_sorted:?}");
        let slp_res = apply_eij_poly_on_program::<_,_,9>(&slp, &casimir_sorted, i64_to_c).unwrap();

        // println!("SLP Res: \n{}", stringify_slp(&slp_res));
        println!("SLP Res: \n{}", stepwise_slp_to_poly(&slp_res, Rational64::ONE).split("\n").last().unwrap());


    }

    #[test]
    fn test_apply_casimir_on_program_reduced() {
        let slp_str = "=C<2,0,0>
=C<0,2,0>
=C<0,0,2>
=L0*L1
=L3*L2";
        let slp = slp_parser_rational(slp_str).expect("Should parse SLP");
        let i64_to_c = |i| Rational64::new(i, 1);
        let casimir_sorted = eval_proj_pairs_to_sorted_basis(&vec![(2,-12),(2,-12),(2,-28)], 3);

        println!("Casimir el is {casimir_sorted:?}");
        let slp_res = apply_eij_poly_on_program::<_,_,9>(&slp, &casimir_sorted, i64_to_c).unwrap();

        // println!("SLP Res: \n{}", stringify_slp(&slp_res));
        println!("SLP Res: \n{}", stepwise_slp_to_poly(&slp_res, Rational64::ONE).split("\n").last().unwrap());
        // println!("SLP Res: \n{}", stepwise_slp_to_poly(&slp_res, Rational64::ONE));
    }

    #[test]
    fn test_example_5_4_projection() {
        let slp_str = "=C<2,0,0>
=C<0,2,0>
=C<0,0,2>
=L0*L1
=L3*L2";
        let slp = slp_parser_rational(slp_str).expect("Should parse SLP");
        let i64_to_c = |i| Rational64::new(i, 1);

        let slp_res = apply_lambda_projection_to_slp::<_,_,9>(&slp, &vec![6,0,0], i64_to_c).unwrap();
        println!("SLP Res: \n{}", stepwise_slp_to_poly(&slp_res, Rational64::ONE).split("\n").last().unwrap());
    }
}

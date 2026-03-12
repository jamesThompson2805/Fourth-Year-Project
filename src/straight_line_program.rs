//! straight_line_program includes the type definition for a straight line program and associated objects and operations.
//! Type definitions are of metavariables, line reference types, expressions (either metavariable or summation or product) and a straight line program.
//! Operations include checks for correctness that cannot be verified by type, operations performed on the SLP and evaluation at a point.
#![allow(dead_code)]

/// MetaVariables refer to a monomial in the underlying polynomial, characterised here by an ordering of the variables and the degree sequence of the monomial.
/// A vector is used as we cannot know the number of variables in the underlying polynomial (often denoted k) at compile time.
pub type MetaVar = Vec<u32>;

/// SLPVar is a variable in a summation or product expression, it can be either a numeric type or a line reference.
///  SLPVar is generic as we don't apply any actual addition or multiplication upon the coefficients
#[derive(Debug, Clone, Copy, Hash, Eq)]
pub enum SLPVar<T> {
    L(usize),
    C(T),
}
impl<T> PartialEq for SLPVar<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (SLPVar::L(n), SLPVar::L(m)) => n == m,
            (SLPVar::C(p), SLPVar::C(q)) => p == q,
            _ => false,
        }
    }
}

/// Operation is either addition or multiplication, used in compound expressions (eg. L1*L3)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Operation {
    Plus,
    Mult,
}

/// SLPLine is either a metavariable or an addition or multiplication expression (tuple of three elements)
#[derive(Debug, Clone, Hash, Eq)]
pub enum SLPLine<T> {
    Input(MetaVar),
    Compound((SLPVar<T>, Operation, SLPVar<T>)),
}
impl<T> PartialEq for SLPLine<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (SLPLine::Input(m1), SLPLine::Input(m2)) => m1 == m2,
            (SLPLine::Compound(c1), SLPLine::Compound(c2)) => c1 == c2,
            _ => false,
        }
    }
}

/// SLP is a vector of lines
pub type SLP<T> = Vec<SLPLine<T>>;

use std::fmt::Display;
/// stringify_slpvar returns a string expression of a slp variable
fn stringify_slpvar<T: Display>(var: &SLPVar<T>) -> String {
    use SLPVar::*;
    match var {
        C(c) => c.to_string(),
        L(lno) => format!("L{lno}"),
    }
}
/// stringify_slpline returns a string expression of a slp line
pub fn stringify_slpline<T: Display>(line: &SLPLine<T>) -> String {
    use SLPLine::*;
    match line {
        Input(m) => format!(
            "C<{}>",
            m.into_iter()
                .map(|i| i.to_string())
                .collect::<Vec<String>>()
                .join(",")
        ),
        Compound((s1, Operation::Plus, s2)) => {
            format!("{} + {}", stringify_slpvar(s1), stringify_slpvar(s2))
        }
        Compound((s1, Operation::Mult, s2)) => {
            format!("{} * {}", stringify_slpvar(s1), stringify_slpvar(s2))
        }
    }
}
/// stringify_slp returns a string expression of a slp
pub fn stringify_slp<T: Display>(p_slp: &SLP<T>) -> String {
    p_slp
        .iter()
        .enumerate()
        .map(|(i, l)| {
            "L".to_string() + i.to_string().as_str() + ": " + stringify_slpline(l).as_str()
        })
        .collect::<Vec<String>>()
        .join("\n")
}

/// sum_slp consumes two straight line programs and returns a new one that evaluates to the sum of the inputs
/// it may duplicate metavariable definition
pub fn add_slp<T>(mut slp1: SLP<T>, mut slp2: SLP<T>) -> SLP<T> {
    use Operation::*;
    use SLPLine::*;
    use SLPVar::*;
    let len1 = slp1.len();
    let len2 = slp2.len();
    let mut slp: Vec<SLPLine<T>> = vec![];
    slp.append(&mut slp1);
    slp.append(&mut slp2);
    for line in slp.iter_mut().skip(len1) {
        match line {
            Compound((L(n1), _, L(n2))) => {
                *n1 += len1;
                *n2 += len1
            }
            Compound((L(n), _, _)) | Compound((_, _, L(n))) => *n += len1,
            _ => (),
        }
    }
    slp.push(Compound((L(len1 - 1), Plus, L(len1 + len2 - 1))));
    slp
}
/// mult_slp consumes two straight line programs and returns a new one that evaluates to the product of the inputs
/// it may duplicate metavariable definition
pub fn mult_slp<T>(mut slp1: SLP<T>, mut slp2: SLP<T>) -> SLP<T> {
    use Operation::*;
    use SLPLine::*;
    use SLPVar::*;
    let len1 = slp1.len();
    let len2 = slp2.len();
    let mut slp: Vec<SLPLine<T>> = vec![];
    slp.append(&mut slp1);
    slp.append(&mut slp2);
    for line in slp.iter_mut().skip(len1) {
        match line {
            Compound((L(n1), _, L(n2))) => {
                *n1 += len1;
                *n2 += len1
            }
            Compound((L(n), _, _)) | Compound((_, _, L(n))) => *n += len1,
            _ => (),
        }
    }
    slp.push(Compound((L(len1 - 1), Mult, L(len1 + len2 - 1))));
    slp
}

/// scale_slp modifies a SLP to evaluates to the original slp scaled by coefficient c
pub fn scale_slp<T>(slp: &mut SLP<T>, c: T) {
    use Operation::*;
    use SLPLine::*;
    use SLPVar::*;
    let len1 = slp.len();
    slp.push(Compound((L(len1 - 1), Mult, C(c))));
}

use std::collections::HashSet;
/// get_metavars_used returns a set of all metavars used in the slp
pub fn get_metavars_used<T>(slp: &SLP<T>) -> HashSet<MetaVar> {
    use SLPLine::*;
    slp.iter()
        .filter_map(|l| match l {
            Input(m) => Some(m.clone()),
            _ => None,
        })
        .collect()
}

use itertools::Itertools;
/// are_metavars_all_same_len returns whether all metavars in slp are same length
/// this is an important validity check for a SLP to pass
pub fn are_metavars_all_same_len<T>(slp: &SLP<T>) -> bool {
    get_metavars_used(slp).iter().map(|m| m.len()).all_equal()
}
/// get_metavar_len returns Some(length) if metavars are present and all have same length, otherwise None.
pub fn get_metavar_len<T>(slp: &SLP<T>) -> Option<usize> {
    match are_metavars_all_same_len(slp) {
        true => get_metavars_used(slp).iter().next().map(|m| m.len()),
        false => None,
    }
}

/// do_metavars_refer_to_homogeneous_poly returns whether all metavars refer to monomials of same degree
/// this is another important validity check for a SLP to pass
pub fn do_metavars_refer_to_homogeneous_poly<T>(slp: &SLP<T>) -> bool {
    get_metavars_used(slp)
        .iter()
        .map(|m| m.into_iter().sum::<u32>())
        .all_equal()
}

/// get_metavar_monomial_degree gets the degree of monomial that all metavars refer to
pub fn get_metavar_monomial_degree<T>(slp: &SLP<T>) -> Option<u32> {
    match do_metavars_refer_to_homogeneous_poly(slp) {
        true => get_metavars_used(slp).iter().next().map(|m| m.iter().sum()),
        false => None,
    }
}

/// are_all_line_refs_valid returns option of whether all references to other lines in the program come before the current line
///  if there are invalid lines, the first is provided
///  this enforces the property that the SLP is a Directed Acyclic Graph
pub fn are_all_line_refs_valid<T>(slp: &SLP<T>) -> Option<usize> {
    use SLPLine::*;
    use SLPVar::*;
    slp.iter()
        .enumerate()
        .find(|(i, line)| match line {
            Compound((L(n1), _, L(n2))) => n1 >= &i || n2 >= &i,
            Compound((L(n), _, _)) | Compound((_, _, L(n))) => n >= &i,
            Compound(_) => false,
            Input(_) => false,
        })
        .map(|(i, _)| i)
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
enum SLPLineReduced<T> {
    L(SLPLine<T>),
    V(SLPVar<T>),
}

use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, Mul};
/// reduce_slp shrinks the size of the current slp whilst reducing none of the information it
/// contains, by removing obvious duplicated lines and recoordinating the line references.
///
/// slp: The slp to be reduced.
/// zero: the value of zero in the custom data type T.
/// one: The value of one in the custom data type T.
pub fn reduce_slp<T>(slp: &mut SLP<T>, zero: T, one: T) -> Vec<usize>
where
    T: Add<Output = T> + Mul<Output = T> + Clone + Eq + Hash + Debug,
{
    let mut similar_lines: HashMap<usize, (usize, SLPLineReduced<T>)> = HashMap::new();
    let mut lines_seen: HashMap<SLPLineReduced<T>, usize> = HashMap::new();

    use Operation::*;
    use SLPLine::*;
    use SLPLineReduced::*;
    use SLPVar::C;
    use SLPVar::L as VL;
    for (lno, line) in slp.iter().enumerate() {
        let slpline_reduced = match line {
            Input(v) => L(Input(v.clone())),
            Compound((C(c1), Plus, C(c2))) => V(C(c1.clone() + c2.clone())),
            Compound((C(c1), Mult, C(c2))) => V(C(c1.clone() * c2.clone())),
            Compound((C(c), Plus, VL(n))) | Compound((VL(n), Plus, C(c))) => {
                let term = similar_lines.get(n).map_or(VL(*n), |p| {
                    if let V(var) = &p.1 {
                        var.clone()
                    } else {
                        VL(p.0)
                    }
                });
                match term {
                    VL(n) => {
                        if c == &zero {
                            V(term)
                        } else {
                            L(Compound((C(c.clone()), Plus, VL(n))))
                        }
                    }
                    C(c2) => V(C(c2 + c.clone())),
                }
            }
            Compound((C(c), Mult, VL(n))) | Compound((VL(n), Mult, C(c))) => {
                let term = similar_lines.get(n).map_or(VL(*n), |p| {
                    if let V(var) = &p.1 {
                        var.clone()
                    } else {
                        VL(p.0)
                    }
                });
                match term {
                    VL(n) => {
                        if c == &zero {
                            V(C(c.clone()))
                        } else if c == &one {
                            V(term)
                        } else {
                            L(Compound((C(c.clone()), Mult, VL(n))))
                        }
                    }
                    C(c2) => V(C(c2 * c.clone())),
                }
            }
            Compound((VL(n1), Plus, VL(n2))) => {
                let term1 = similar_lines.get(n1).map_or(VL(*n1), |p| {
                    if let V(var) = &p.1 {
                        var.clone()
                    } else {
                        VL(p.0)
                    }
                });
                let term2 = similar_lines.get(n2).map_or(VL(*n2), |p| {
                    if let V(var) = &p.1 {
                        var.clone()
                    } else {
                        VL(p.0)
                    }
                });
                match (term1, term2) {
                    (C(c1), C(c2)) => V(C(c1 + c2)),
                    (C(c), VL(n)) | (VL(n), C(c)) => {
                        if c == zero {
                            V(VL(n))
                        }
                        // covered case that coefficient is zero again
                        else {
                            L(Compound((C(c), Plus, VL(n))))
                        }
                    }
                    (VL(n1), VL(n2)) => L(Compound((VL(n1.min(n2)), Plus, VL(n1.max(n2))))),
                }
            }
            Compound((VL(n1), Mult, VL(n2))) => {
                let term1 = similar_lines.get(n1).map_or(VL(*n1), |p| {
                    if let V(var) = &p.1 {
                        var.clone()
                    } else {
                        VL(p.0)
                    }
                });
                let term2 = similar_lines.get(n2).map_or(VL(*n2), |p| {
                    if let V(var) = &p.1 {
                        var.clone()
                    } else {
                        VL(p.0)
                    }
                });
                match (term1, term2) {
                    (C(c1), C(c2)) => V(C(c1 * c2)),
                    (C(c), VL(n)) | (VL(n), C(c)) => {
                        if c == zero {
                            V(C(zero.clone()))
                        }
                        // covered case that coefficient is zero
                        else if c == one {
                            V(VL(n))
                        }
                        // covered case that coefficient is one
                        else {
                            L(Compound((C(c), Mult, VL(n))))
                        }
                    }
                    (VL(n1), VL(n2)) => L(Compound((VL(n1.min(n2)), Mult, VL(n1.max(n2))))),
                }
            }
        };
        // println!("Line {lno} becomes {slpline_reduced:?}");
        if let Some(n) = lines_seen.get(&slpline_reduced) {
            // println!("  Line {lno} reduces to {n}");
            similar_lines.insert(lno, (*n, slpline_reduced));
        } else {
            similar_lines.insert(lno, (lno, slpline_reduced.clone()));
            lines_seen.insert(slpline_reduced, lno);
            lines_seen.insert(V(VL(lno)), lno);
        }
    }

    // find those lines to remove (those that are similar to an earlier line)
    let mut lines_to_remove = similar_lines
        .iter()
        .filter(|(k, v)| **k != v.0)
        .map(|(k, _)| k)
        .collect::<HashSet<_>>();
    lines_to_remove.remove(&(slp.len() - 1)); // should ensure we keep the final line

    // the index_map maps line references to their position after deleting all lines in lines_to_remove
    let mut index_map = vec![0; slp.len()];
    let mut i = 0;
    for j in 0..slp.len() {
        while lines_to_remove.contains(&i) {
            i += 1;
        }
        if i >= slp.len() {
            break;
        }
        index_map[i] = j;
        i += 1;
    }

    // update all references in the slp
    for line in slp.iter_mut() {
        match line {
            Compound((VL(n1), _, VL(n2))) => {
                *n1 = index_map[similar_lines[n1].0];
                *n2 = index_map[similar_lines[n2].0];
            }
            Compound((VL(n), _, _)) | Compound((_, _, VL(n))) => *n = index_map[similar_lines[n].0],
            _ => (),
        }
    }

    let mut current_idx = 0;

    slp.retain_mut(|_line| {
        let should_keep = !lines_to_remove.contains(&current_idx);

        current_idx += 1;
        should_keep
    });

    for (k, v) in similar_lines.iter() {
        index_map[*k] = index_map[v.0];
    }
    // println!("  index map {index_map:?}");
    index_map
}

fn generate_all_d_k_metavars(d: usize, k: usize) -> Vec<Vec<u32>> {
    if k == 0 {
        return vec![];
    } else if k == 1 {
        return vec![vec![d as u32]];
    } else if d == 0 {
        return vec![vec![0;k]];
    }
    (0..=d)
        .map(|i| (i, generate_all_d_k_metavars(d - i, k - 1)))
        .map(|(i, mut vs)| {
            vs.iter_mut().for_each(|v| v.push(i as u32));
            vs
        })
        .flatten()
        .collect()
}

use rand::seq::IteratorRandom;

/// generate_random_slp generates a new random slp using the parameters given.
///
/// num_gates is the number of gates used in the slp
/// gen_coeff is a mutating function that generates a random coefficient
/// d is the degree of the underlying polynomial the metapolynomial references
/// k is the number of variables in the underlying polynomial
/// prob_choose_mult is a 0-1 float representing the probability of a compound gate being Multiply instead of add
/// prob_choose_ref is a 0-1 float representing the probability of a compound gate using a line reference instead of a float
/// prob_choose_metavar is a 0-1 float of prob choosing metavar over ref to compound line
pub fn generate_random_slp<T, F, R>(
    num_gates: usize,
    mut gen_coeff: F,
    d: usize,
    k: usize,
    prob_choose_mult: f64,
    prob_choose_ref: f64,
    prob_choose_metavar: f64,
    rng: &mut R,
) -> SLP<T>
where
    F: FnMut() -> T,
    R: rand::Rng,
{
    use Operation::*;
    use SLPLine::*;
    use SLPVar::*;
    // Steps:
    // 1. generate slp with all valid metavariables already present
    // 2. whilst not added num_gates to slp:
    //      2.1 choose gate type (using prob_mult)
    //      2.2 choose the variables (using prob_choose_ref)
    //      2.3 choose the line reference (using prob_choose_metavar)
    // 1. Generate slp with all valid metavariables already present
    let all_metavars = generate_all_d_k_metavars(d, k);
    let all_metavars_len = all_metavars.len();
    let mut unseen_lines: HashSet<usize> = HashSet::new();
    let mut slp: SLP<T> = Vec::new();
    for m in all_metavars {
        slp.push(Input(m));
    }

    for _ in 0..num_gates {
        let term1 = if rng.random::<f64>() < prob_choose_ref {
            if slp.len() > all_metavars_len && rng.random::<f64>() > prob_choose_metavar {
                let line_to_remove = *unseen_lines.iter().choose(rng).unwrap_or( &rng.random_range(all_metavars_len..slp.len()) );
                unseen_lines.remove( &line_to_remove );
                L( line_to_remove )
            } else {
                L(rng.random_range(0..all_metavars_len))
            }
        } else {
            C(gen_coeff())
        };
        let term2 = if rng.random::<f64>() < prob_choose_ref {
            if slp.len() > all_metavars_len && rng.random::<f64>() > prob_choose_metavar {
                let line_to_remove = *unseen_lines.iter().choose(rng).unwrap_or( &rng.random_range(all_metavars_len..slp.len()) );
                unseen_lines.remove( &line_to_remove );
                L( line_to_remove )
            } else {
                L(rng.random_range(0..all_metavars_len))
            }
        } else {
            C(gen_coeff())
        };
        if rng.random::<f64>() < prob_choose_mult {
            slp.push(Compound((term1, Mult, term2)));
        } else {
            slp.push(Compound((term1, Plus, term2)));
        }
        unseen_lines.insert(slp.len()-1);
    }

    slp
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Rational64;
    use rand::{Rng, rng};
    use crate::evaluation::stepwise_slp_to_poly;

    #[test]
    fn gen_random_slp() {
        let mut coeff_rng = rng();

        let gen_coeff = || {
            let numer: i64 = coeff_rng.random_range(-100..100);
            let denom: i64 = coeff_rng.random_range(1..100);  // strictly positive to avoid division by zero
            Rational64::new(numer, denom)
        };

        let num_gates = 10;
        let d = 3;
        let k = 3;
        let slp = generate_random_slp::<Rational64,_,_>(num_gates, gen_coeff, d, k, 0.6, 0.8, 0.3, &mut rng());
        // println!("{}", stepwise_slp_to_poly(&slp, Rational64::ONE).split("\n").last().unwrap());
        println!("{}", stepwise_slp_to_poly(&slp, Rational64::ONE));
    }

}
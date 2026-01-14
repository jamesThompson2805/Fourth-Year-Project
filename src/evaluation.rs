//! evaluation.rs is dedicated to functions assessing if SLPs are the same or not
#![allow(dead_code)]

use super::straight_line_program::*;

use rand::Rng;
use rand::distr::Distribution;

use std::collections::BTreeMap;
use std::collections::HashMap;
use std::ops::{Add,Mul};
use std::fmt::Display;

use num_rational::Rational64;
use num_complex::Complex64;
/// get_random_val_for_metavars is a function that generates a random point to evaluate a SLP at by generating a coefficient value for each metavariable
/// T: The coefficient type
/// D: The Generator object that takes the rng and converts it into a coefficient
/// R: The source of randomness
/// returns a map of metavariables to a point representing its value
pub fn get_random_val_for_metavars<T, D, R>(slp: &SLP<T>, generator: &D, rng: &mut R) -> HashMap<MetaVar,T>
where 
    D: Distribution<T>,
    R: Rng + ?Sized,
{
    // The generator "samples" the randomness to produce T
    let mut metavar_points: HashMap<MetaVar,T> = HashMap::new();
    for m in get_metavars_used(slp) {
        metavar_points.insert(m, generator.sample(rng));
    }
    metavar_points
}

/// eval_slp_at_point takes a dictionary of metavariables and choice of value pairs and a SLP
/// returns the evaluated SLP if the dictionary contains all necessary metavariables otherwise None
pub fn eval_slp_at_point<T: Add<Output=T> + Mul<Output=T> + Copy>(input_vals: &HashMap<MetaVar, T>, slp: &SLP<T>) -> Option<T> {
    use SLPLine::*;
    use SLPVar::*;
    use Operation::*;
    let mut line_evals: Vec<T> = vec![];
    for line in slp {
        match line {
            Input(m) => line_evals.push( *(input_vals.get(m)?) ),
            Compound((s1, op, s2)) => {
                let s1_val: T = *match s1 { C(c) => c, L(n) => &line_evals[*n]};
                let s2_val: T = *match s2 { C(c) => c, L(n) => &line_evals[*n]};

                match op {
                    Plus => line_evals.push( s1_val + s2_val),
                    Mult => line_evals.push( s1_val * s2_val),
                }
            },
        }
    }
    line_evals.last().copied()
}

/// are_slps_similar returns if two slps match at for a given selection of input values
pub fn are_slps_similar<T: Add<Output=T> + Mul<Output=T> + PartialEq + Copy>(input_vals: &HashMap<MetaVar, T>, slp1: &SLP<T>, slp2: &SLP<T>) -> bool {
    eval_slp_at_point(input_vals, slp1) == eval_slp_at_point(input_vals, slp2)
}

/// RationalDistr is a distribution for generating random rational numbers
pub struct RationalDistr { }
impl Distribution<Rational64> for RationalDistr {
    // generates rational number in range [-10,10]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Rational64 {
        Rational64::new( rng.random_range(-10..=10), rng.random_range(1..=10) )
    }
}

/// ComplexDistr is a distribution for generating random complex numbers
pub struct ComplexDistr { }
impl Distribution<Complex64> for ComplexDistr {
    // generates complex number with standard uniform components
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Complex64 {
        Complex64::new( rng.random(), rng.random() )
    }
}

/// are_complex_slp_similar takes two SLP with complex coefficients and evaluates if similar on random point
pub fn are_complex_slp_similar(slp1: &SLP<Complex64>, slp2: &SLP<Complex64>) -> (HashMap<MetaVar,Complex64>, bool) {
    let metavars1 = get_metavars_used(slp1);
    let metavars2 = get_metavars_used(slp2);
    let metavars = metavars1.union( &metavars2 );
    let mut rng = rand::rng();
    let dist = ComplexDistr {};

    let input_vals = metavars.map(|m| (m.clone(), dist.sample(&mut rng) )).collect();
    let similar = are_slps_similar(&input_vals, slp1, slp2);
    (input_vals, similar)
}
/// are_real_slp_similar takes two SLP with f64 coefficients and evaluates if similar on random point
pub fn are_real_slp_similar(slp1: &SLP<f64>, slp2: &SLP<f64>) -> (HashMap<MetaVar,f64>, bool) {
    let metavars1 = get_metavars_used(slp1);
    let metavars2 = get_metavars_used(slp2);
    let metavars = metavars1.union( &metavars2 );
    let mut rng = rand::rng();

    let input_vals = metavars.map(|m| (m.clone(), rng.random_range(-10.0..=10.0)) ).collect();
    let similar = are_slps_similar(&input_vals, slp1, slp2);
    (input_vals, similar)
}
/// are_rational_slp_similar takes two SLP with rational coefficients and evaluates if similar on random point
pub fn are_rational_slp_similar(slp1: &SLP<Rational64>, slp2: &SLP<Rational64>) -> (HashMap<MetaVar,Rational64>, bool) {
    let metavars1 = get_metavars_used(slp1);
    let metavars2 = get_metavars_used(slp2);
    let metavars = metavars1.union( &metavars2 );
    let mut rng = rand::rng();
    let dist = RationalDistr {};

    let input_vals = metavars.map(|m| (m.clone(), dist.sample(&mut rng) )).collect();
    let similar = are_slps_similar(&input_vals, slp1, slp2);
    (input_vals, similar)
}

type Poly<T> = BTreeMap< BTreeMap<Vec<u32>, u32>, T>;

fn add_poly<T> (p1: Poly<T>, p2: &Poly<T>) -> Poly<T>
where T: Add<Output=T> + Default + PartialEq + Clone,
{
    let mut res = p1.clone();
    for term in p2 {
        // term in both polynomials => add coefficients
        if let Some(x) = res.get_mut(term.0) {
            if x.clone() + term.1.clone() == T::default() {
                res.remove(&term.0);
            } else {
                *x = x.clone() + term.1.clone();
            }
        } else { // term in second poly only => add term alongside coefficient
            if !(term.1.clone() == T::default()) {
                res.insert(term.0.clone(), term.1.clone());
            }
        }
    }
    res
}

fn mul_mono<T>(p1: &Poly<T>, m2: (&BTreeMap<Vec<u32>, u32>,&T)) -> Poly<T>
where T: Mul<Output=T> + Default + PartialEq + Clone,
{
    if m2.0.is_empty() || m2.1 == &T::default() {
        return BTreeMap::new();
    }

    let mut res = BTreeMap::new();
    for term in p1 {
        let mut res_term = term.0.clone();

        for var in m2.0 {
            if let Some(x) = res_term.get_mut(var.0) {
                *x = x.clone() + *var.1;
            } else {
                res_term.insert(var.0.clone(), *var.1);
            }
        }
        res.insert(res_term, term.1.clone() * m2.1.clone());
    }

    res
}

fn mul_poly<T>(p1: &Poly<T>, p2: &Poly<T>) -> Poly<T>
where T: Add<Output=T> + Mul<Output=T> + Default + PartialEq + Clone,
{
    let terms: Vec<Poly<T>> = p2.iter().map(|t| mul_mono(&p1, t)).collect();
    if terms.len() == 0 {
        BTreeMap::new()
    } else {
        let first = terms[0].clone();
        terms.iter().skip(1).fold(first, |acc, el| add_poly(acc, el))
    }
}

fn scale_poly<T>(p1: &Poly<T>, c: T) -> Poly<T>
where T: Mul<Output=T> + Clone + PartialEq + Default,
{
    if c == T::default() {
        return BTreeMap::new();
    }
    p1.iter().map(|(k,v)| (k.clone(), v.clone() * c.clone() )).collect()
}

fn stringify_pvar(v: &Vec<u32>) -> String {
    "C<".to_string() + &v.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",") + ">"
}
fn stringify_mono(m: &BTreeMap<Vec<u32>,u32>) -> String {
    m.iter().map(|t| stringify_pvar(t.0) + "^" + &t.1.to_string()).collect::<Vec<_>>().join(".")
}
fn stringify_poly<T: Display>(p: &Poly<T>) -> String {
    p.iter().map(|(m,c)| "(".to_string() + &c.to_string() + ")" + &stringify_mono(m)).collect::<Vec<_>>().join(" + ")
}

pub fn stepwise_slp_to_poly<T>(slp: &SLP<T>, unit: T) -> String
where T: Add<Output=T> + Mul<Output=T> + Default + PartialEq + Clone + Display,
{
    use SLPLine::*;
    use SLPVar::*;
    use Operation::*;
    let mut lines: Vec<Poly<T>> = Vec::new();
    for line in slp {
        lines.push( 
        match line {
            Input(m) => BTreeMap::from([(  BTreeMap::from([(m.clone(), 1)]), unit.clone()  )]),

            Compound(( C(t1), Plus, C(t2) )) => {
                if t1.clone() + t2.clone() == T::default() { // handle case of value being 0
                    BTreeMap::new()
                } else {
                    BTreeMap::from([( BTreeMap::new(), t1.clone() + t2.clone() )])
                }
            },
            Compound(( C(t1), Mult, C(t2) )) => {
                if t1.clone() * t2.clone() == T::default() { // handle case of value being 0
                    BTreeMap::new()
                } else {
                    BTreeMap::from([( BTreeMap::new(), t1.clone() * t2.clone() )])
                }
            },

            Compound(( C(t), Plus, L(n) )) | Compound(( L(n), Plus, C(t) )) => {
                let mut res = lines[*n].clone();
                if let Some(x) = res.get_mut(&BTreeMap::new() ) {
                    *x = x.clone() + t.clone();
                } else {
                    res.insert( BTreeMap::new(), t.clone() );
                }
                res
            },
            Compound(( C(t), Mult, L(n) )) | Compound(( L(n), Mult, C(t) )) => {
                scale_poly(&lines[*n], t.clone())
            },
            Compound(( L(n1), Plus, L(n2) )) => add_poly( lines[*n1].clone(), &lines[*n2] ),
            Compound(( L(n1), Mult, L(n2) )) => mul_poly( &lines[*n1],&lines[*n2] ),
        });
    }

    (0..slp.len()).map(
        |i| format!("L{i: <4}: {: <15} : {}", stringify_slpline(&slp[i]), stringify_poly(&lines[i]) )
    ).collect::<Vec<_>>().join("\n")
}
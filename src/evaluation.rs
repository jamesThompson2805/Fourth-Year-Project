//! evaluation.rs is dedicated to functions assessing if SLPs are the same or not
#![allow(dead_code)]

use super::straight_line_program::*;

use rand::Rng;
use rand::distr::Distribution;
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

use std::collections::HashMap;
use std::ops::Add;
use std::ops::Mul;
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

use num_rational::Rational64;
/// RationalDistr is a distribution for generating random rational numbers
pub struct RationalDistr { }
impl Distribution<Rational64> for RationalDistr {
    // generates rational number in range [-10,10]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Rational64 {
        Rational64::new( rng.random_range(-10..=10), rng.random_range(1..=10) )
    }
}

use num_complex::Complex64;
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
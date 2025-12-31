//! straight_line_program includes the type definition for a straight line program and associated objects and operations.
//! Type definitions are of metavariables, line reference types, expressions (either metavariable or summation or product) and a straight line program.
//! Operations include checks for correctness that cannot be verified by type, operations performed on the SLP and evaluation at a point.
#![allow(dead_code)]

/// MetaVariables refer to a monomial in the underlying polynomial, characterised here by an ordering of the variables and the degree sequence of the monomial.
/// A vector is used as we cannot know the number of variables in the underlying polynomial (often denoted k) at compile time. 
pub type MetaVar = Vec<u32>;

/// SLPVar is a variable in a summation or product expression, it can be either a numeric type or a line reference.
///  SLPVar is generic as we don't apply any actual addition or multiplication upon the coefficients
#[derive(Debug, Clone, Copy)]
pub enum SLPVar<T> { L(usize), C(T) }


/// Operation is either addition or multiplication, used in compound expressions (eg. L1*L3)
#[derive(Debug, Clone, Copy)]
pub enum Operation { Plus, Mult, }

/// SLPLine is either a metavariable or an addition or multiplication expression (tuple of three elements)
#[derive(Debug, Clone)]
pub enum SLPLine<T> {
    Input(MetaVar),
    Compound((SLPVar<T>, Operation, SLPVar<T>)),
}

/// SLP is a vector of lines
pub type SLP<T> = Vec<SLPLine<T>>;

use std::fmt::Display;
/// stringify_slpvar returns a string expression of a slp variable
fn stringify_slpvar<T:Display>(var: &SLPVar<T>) -> String {
    use SLPVar::*;
    match var {
        C(c) => c.to_string(),
        L(lno) => format!("L{lno}"),
    }
}
/// stringify_slpline returns a string expression of a slp line
fn stringify_slpline<T: Display>(line: &SLPLine<T>) -> String {
    use SLPLine::*;
    match line {
        Input(m) => format!("C<{}>",m.into_iter().map(|i| i.to_string()).collect::<Vec<String>>().join(",")),
        Compound((s1,Operation::Plus, s2)) => format!("{} + {}", stringify_slpvar(s1), stringify_slpvar(s2)),
        Compound((s1,Operation::Mult, s2)) => format!("{} * {}", stringify_slpvar(s1), stringify_slpvar(s2)),
    }
}
/// stringify_slp returns a string expression of a slp
pub fn stringify_slp<T: Display>(p_slp: &SLP<T>) -> String {
    p_slp.iter().enumerate().map(|(i,l)| "L".to_string() + i.to_string().as_str() + ": " + stringify_slpline(l).as_str() )
        .collect::<Vec<String>>().join("\n")
}

/// sum_slp consumes two straight line programs and returns a new one that evaluates to the sum of the inputs
/// it may duplicate metavariable definition
pub fn add_slp<T>(mut slp1: SLP<T>, mut slp2: SLP<T>) -> SLP<T> {
    use SLPLine::*;
    use SLPVar::*;
    use Operation::*;
    let len1 = slp1.len();
    let len2 = slp2.len();
    let mut slp: Vec<SLPLine<T>> = vec![];
    slp.append(&mut slp1);
    slp.append(&mut slp2);
    slp.push(
        Compound(( L(len1 - 1), Plus, L(len1 + len2 - 1)))
    );
    slp
}
/// mult_slp consumes two straight line programs and returns a new one that evaluates to the product of the inputs
/// it may duplicate metavariable definition
pub fn mult_slp<T>(mut slp1: SLP<T>, mut slp2: SLP<T>) -> SLP<T> {
    use SLPLine::*;
    use SLPVar::*;
    use Operation::*;
    let len1 = slp1.len();
    let len2 = slp2.len();
    let mut slp: Vec<SLPLine<T>> = vec![];
    slp.append(&mut slp1);
    slp.append(&mut slp2);
    slp.push(
        Compound(( L(len1 - 1), Mult, L(len1 + len2 - 1)))
    );
    slp
}

/// scale_slp modifies a SLP to evaluates to the original slp scaled by coefficient c
pub fn scale_slp<T>(slp: &mut SLP<T>, c: T) {
    use SLPLine::*;
    use SLPVar::*;
    use Operation::*;
    let len1 = slp.len();
    slp.push(
        Compound(( L(len1 - 1), Mult, C(c)))
    );
}

use std::collections::HashSet;
/// get_metavars_used returns a set of all metavars used in the slp
pub fn get_metavars_used<T>(slp: &SLP<T>) -> HashSet<MetaVar> {
    use SLPLine::*;
    slp.iter().filter_map(|l| match l { Input(m) => Some(m.clone()), _=>None }).collect()
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
    get_metavars_used(slp).iter().map(|m| m.into_iter().sum::<u32>()).all_equal()
}

/// get_metavar_monomial_degree gets the degree of monomial that all metavars refer to
pub fn get_metavar_monomial_degree<T>(slp: &SLP<T>) -> Option<u32> {
    match do_metavars_refer_to_homogeneous_poly(slp) {
        true => get_metavars_used(slp).iter().next().map(|m| m.iter().sum()),
        false => None,
    }
}

/// are_all_line_refs_valid returns boolean of whether all references to other lines in the program come before the current line
///  this enforces the property that the SLP is a Directed Acyclic Graph
pub fn are_all_line_refs_valid<T>(slp: &SLP<T>) -> bool {
    use SLPLine::*;
    use SLPVar::*;
    slp.iter().enumerate().all(|(i,line)| match line {
        Compound(( L(n1), _, L(n2) )) => n1 < &i && n2 < &i,
        Compound(( L(n), _, _ )) | Compound(( _, _, L(n) )) => n < &i,
        Compound(_) => true,
        Input(_) => true,
    })
}


/* 
use std::collections::BTreeMap;
type MetaMonomial = BTreeMap<MetaVar, u64>;
type MetaPolynomial = BTreeMap<MetaMonomial, Complex64>;

fn create_coeff(c: &Complex64) -> MetaPolynomial {
    BTreeMap::from([( BTreeMap::from([(Vec::new(),1)]), *c )])
}

fn evaluate_plus( s1: &SLPVar, s2: &SLPVar, so_far: &Vec<Option<MetaPolynomial>>) -> Option<MetaPolynomial> {
    use SLPVar::*;
    let mp_1 = match s1 {
        L(n) => so_far[*n].clone(),
        F(c) => create_coeff(c),
    };
    let mut mp_2 = match s2 {
        L(n) => so_far[*n].clone(),
        F(c) => create_coeff(c),
    };
    let mut mp: MetaPolynomial = BTreeMap::new();
    for (k,v) in mp_1 {
        if mp_2.contains_key(&k) {
            let v_new = v + mp_2[&k];
            mp_2.remove(&k);
            mp.insert(k, v_new);
        } else {
            mp.insert(k,v);
        }
    }
    for (k,v) in mp_2 {
        mp.insert(k,v);
    }
    unimplemented!()
}

pub fn evaluate_slp(slp: &SLP) -> String {
    use SLPLine::*;
    use Operation::*;
    let mut mp_v: Vec<Option<MetaPolynomial>> = Vec::new();
    for line in slp {
        let mp = match line {
                Input(m) => {
                    let mm:MetaMonomial = BTreeMap::from([(m.clone(),1)]);
                    Some(BTreeMap::from( [(mm, Complex64::ONE)] ))
                },
                Compound((s1,Plus, s2)) => ,
                Compound((s1,Mult, s2)) => ,

        };
        mp_v.push(mp);
        
    }
    unimplemented!()
}
*/
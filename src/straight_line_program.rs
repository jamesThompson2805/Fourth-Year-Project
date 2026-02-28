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
pub enum SLPVar<T> { L(usize), C(T) }
impl<T> PartialEq for SLPVar<T> where T: PartialEq {
    fn eq(&self, other: &Self) -> bool {
        match (self,other) {
            (SLPVar::L(n),SLPVar::L(m)) => n==m,
            (SLPVar::C(p),SLPVar::C(q)) => p==q,
            _ => false,
        }
    }
}


/// Operation is either addition or multiplication, used in compound expressions (eg. L1*L3)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Operation { Plus, Mult, }

/// SLPLine is either a metavariable or an addition or multiplication expression (tuple of three elements)
#[derive(Debug, Clone, Hash, Eq)]
pub enum SLPLine<T> {
    Input(MetaVar),
    Compound((SLPVar<T>, Operation, SLPVar<T>)),
}
impl<T> PartialEq for SLPLine<T> where T: PartialEq {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (SLPLine::Input(m1), SLPLine::Input(m2)) => m1==m2,
            (SLPLine::Compound(c1), SLPLine::Compound(c2)) => c1==c2,
            _ => false,
        }
    }
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
pub fn stringify_slpline<T: Display>(line: &SLPLine<T>) -> String {
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
    for line in slp.iter_mut().skip(len1) {
        match line {
            Compound(( L(n1), _, L(n2) )) => {*n1 += len1; *n2 += len1},
            Compound(( L(n), _, _ )) | Compound(( _, _, L(n) )) => *n += len1,
            _ => (),
        }
    }
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
    for line in slp.iter_mut().skip(len1) {
        match line {
            Compound(( L(n1), _, L(n2) )) => {*n1 += len1; *n2 += len1},
            Compound(( L(n), _, _ )) | Compound(( _, _, L(n) )) => *n += len1,
            _ => (),
        }
    }
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

/// are_all_line_refs_valid returns option of whether all references to other lines in the program come before the current line
///  if there are invalid lines, the first is provided
///  this enforces the property that the SLP is a Directed Acyclic Graph
pub fn are_all_line_refs_valid<T>(slp: &SLP<T>) -> Option<usize> {
    use SLPLine::*;
    use SLPVar::*;
    slp.iter().enumerate().find(|(i,line)| match line {
        Compound(( L(n1), _, L(n2) )) => n1 >= &i || n2 >= &i,
        Compound(( L(n), _, _ )) | Compound(( _, _, L(n) )) => n >= &i,
        Compound(_) => false,
        Input(_) => false,
    }).map(|(i,_)| i)
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
enum SLPLineReduced<T> {L(SLPLine<T>), V(SLPVar<T>)}

use std::collections::HashMap;
use std::ops::{Add,Mul};
use std::hash::Hash;
pub fn reduce_slp<T>(slp: &mut SLP<T>, zero: T) -> HashMap<usize, usize>
where
    T: Default + Add<Output=T> + Mul<Output=T> + Clone + Eq + Hash,
{
    let mut similar_lines: HashMap<usize, (usize, SLPLineReduced<T>)> = HashMap::new();
    let mut lines_seen: HashMap<SLPLineReduced<T>, usize> = HashMap::new();

    use SLPLine::*;
    use SLPVar::C;
    use SLPVar::L as VL;
    use SLPLineReduced::*;
    use Operation::*;
    for (lno, line) in slp.iter().enumerate() {
        let slpline_reduced = match line {
            Input(v) => L(Input(v.clone())),
            Compound(( C(c1), Plus, C(c2) )) => V(C(c1.clone() + c2.clone())),
            Compound(( C(c1), Mult, C(c2) )) => V(C(c1.clone() * c2.clone())),
            Compound(( C(c), Plus, VL(n) )) | Compound(( VL(n), Plus, C(c) )) => {
                if c == &zero { // optimisation: if plus zero then use only the reference in future
                    let term = similar_lines.get(n).map_or(VL(*n), |p| { if let V(var) = &p.1 {var.clone()} else {VL(p.0)} });
                    V(term)
                } else {
                    L( Compound(( C(c.clone()), Plus, VL(*n) )) )
                }
            },
            Compound(( C(c), Mult, VL(n) )) | Compound(( VL(n), Mult, C(c) )) => {
                if c == &zero { // optimisation: if times zero then map straight to zero reference
                    V(C(c.clone()))
                } else if c == &T::default() { // optimisation: if times one then only use reference in future
                    let term = similar_lines.get(n).map_or(VL(*n), |p| { if let V(var) = &p.1 {var.clone()} else {VL(p.0)} });
                    V(term)
                } else {
                    L( Compound(( C(c.clone()), Mult, VL(*n) )) )
                }
            },
            Compound(( VL(n1), Plus, VL(n2) )) => {
                let term1 = similar_lines.get(n1).map_or(VL(*n1), |p| { if let V(var) = &p.1 {var.clone()} else {VL(p.0)} });
                let term2 = similar_lines.get(n2).map_or(VL(*n2), |p| { if let V(var) = &p.1 {var.clone()} else {VL(p.0)} });
                match (term1, term2) {
                    ( C(c1), C(c2) ) => V(C(c1+c2)),
                    ( C(c), VL(n) ) | ( VL(n), C(c) ) =>  
                        if c == zero { V(VL(n)) } // covered case that coefficient is zero again
                        else { L(Compound(( C(c), Plus, VL(n) ))) }
                    ,
                    ( VL(n1), VL(n2) ) => L(Compound(( VL(n1.min(n2)), Plus, VL(n1.max(n2)) ))),
                }

            },
            Compound(( VL(n1), Mult, VL(n2) )) => {
                let term1 = similar_lines.get(n1).map_or(VL(*n1), |p| { if let V(var) = &p.1 {var.clone()} else {VL(p.0)} });
                let term2 = similar_lines.get(n2).map_or(VL(*n2), |p| { if let V(var) = &p.1 {var.clone()} else {VL(p.0)} });
                match (term1, term2) {
                    ( C(c1), C(c2) ) => V(C(c1*c2)),
                    ( C(c), VL(n) ) | ( VL(n), C(c) ) =>  
                        if c == zero { V(C(zero.clone())) } // covered case that coefficient is zero
                        else if c==T::default() { V(VL(n)) } // covered case that coefficient is one
                        else { L(Compound(( C(c), Mult, VL(n) ))) }
                    ,
                    ( VL(n1), VL(n2) ) => L(Compound(( VL(n1.min(n2)), Mult, VL(n1.max(n2)) ))),
                }
            },
        };
        if let Some(n) = lines_seen.get(&slpline_reduced) {
            similar_lines.insert( lno, (*n, slpline_reduced) );
        } else {
            similar_lines.insert( lno, (lno, slpline_reduced.clone()));
            lines_seen.insert( slpline_reduced, lno);
        }
    }

    // find those lines to remove (those that are similar to an earlier line)
    let lines_to_remove = similar_lines.iter().filter(|(k,v)| **k!=v.0).map(|(k,_)| k).collect::<HashSet<_>>();
    // the index_map maps line references to their position after deleting all lines in lines_to_remove
    let mut index_map = vec![0; slp.len()];
    let mut i=0;
    for j in 0..slp.len() {
        while lines_to_remove.contains(&i) {
            i+=1;
        }
        if i>=slp.len() { break; }
        index_map[i] = j;
        i+=1;
    }

    // update all references in the slp
    for line in slp.iter_mut() {
        match line {
            Compound(( VL(n1), _, VL(n2) )) => {
                *n1 = index_map[ similar_lines[n1].0 ];
                *n2 = index_map[ similar_lines[n1].0 ];
            },
            Compound(( VL(n), _, _ )) | Compound(( _, _, VL(n) )) => *n = index_map[ similar_lines[n].0 ],
            _ => (),
        }
    }

    let mut current_idx = 0;

    slp.retain_mut(|_line| {
        let should_keep = !lines_to_remove.contains(&current_idx);
        
        current_idx += 1;
        should_keep
    });


    similar_lines.into_iter().map(|(k,v)| (k, index_map[v.0])).collect()
}
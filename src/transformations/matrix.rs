//! matrix_transformations provides the singular operation of applying a linear combination of lists of sorted eijs to a slp.
//! This task is considered large enough to dedicate a file to ensure clarity.

#![allow(dead_code)]

use crate::{straight_line_program::*, transformations::apply_eij_on_metavar};

use radix_trie::Trie;
use std::collections::HashMap; 
use radix_trie::TrieKey;
use std::fmt::Display;

/// LinesToTransform holds dict of lines in original SLP and whether they are needed in transformed SLP, with or without transformation
type LinesToTransform<const K: usize> = HashMap<usize, [u32;K]>;
/// LineMap holds dict of lines and transformation applied from original SLP and location in new transformed partial SLP
type LineMap<const K: usize> = HashMap<usize, HashMap<[u32;K],usize>>;
fn insert_lm<const K: usize>(lm: &mut LineMap<K>, ilno: usize, transform: [u32;K], olno: usize) {
    if let Some(hm) = lm.get_mut(&ilno) {
        hm.insert(transform, olno);
    } else {
        lm.insert(ilno, HashMap::from([(transform, olno)]));
    }
}

/// PartialSLPVar is a structure mimicking SLPVar but in an unfinished format allowing it to be converted to a valid SLPVar later
enum PartialSLPVar<T,const K:usize> {
    LIP(usize), // LineInPartial is index of an earlier line in the transformed partial slp
    LTT((usize, [u32;K])), // LineToTranslate, bool represents whether it should be transformed by eij or not
    C(T),
    Zero,
}
/// PartialSLPLine is a structure mimicking SLPLine but in an unfinished format allowing it to be converted to it later
enum PartialSLPLine<T, const K: usize> {
    Input(MetaVar),
    Compound((PartialSLPVar<T,K>, Operation, PartialSLPVar<T,K>)),
}
/// PartialSLPis a structure mimicking SLPbut in an unfinished format allowing it to be converted to it later
type PartialSLP<T, const K: usize> = Vec<PartialSLPLine<T,K>>;

/// stringify_partialslpvar returns string presentation of partialslpvar for debugging
fn stringify_partialslpvar<T: Display, const K: usize>(var: &PartialSLPVar<T,K>) -> String {
    use PartialSLPVar::*;
    match var {
        LTT((lno, transform)) => format!("LTT({lno},T:{transform:?})"),
        C(c) => c.to_string(),
        LIP(lno) => format!("LIP{lno}"),
        Zero => format!("0"),
    }
}
/// stringify_partialslpline returns string presentation of partialslpline for debugging
fn stringify_partialslpline<T: Display, const K: usize>(line: &PartialSLPLine<T,K>) -> String {
    use PartialSLPLine::*;
    match line {
        Input(m) => format!(
            "C<{}>",
            m.into_iter()
                .map(|i| i.to_string())
                .collect::<Vec<String>>()
                .join(",")
        ),
        Compound((s1, Operation::Plus, s2)) => format!(
            "{} + {}",
            stringify_partialslpvar(s1),
            stringify_partialslpvar(s2)
        ),
        Compound((s1, Operation::Mult, s2)) => format!(
            "{} * {}",
            stringify_partialslpvar(s1),
            stringify_partialslpvar(s2)
        ),
    }
}
/// stringify_partialslp returns string presentation of partialslp for debugging
fn stringify_partialslp<T: Display, const K:usize>(p_slp: &PartialSLP<T,K>) -> String {
    p_slp
        .iter()
        .enumerate()
        .map(|(i, l)| {
            "L".to_string() + i.to_string().as_str() + ": " + stringify_partialslpline(l).as_str()
        })
        .collect::<Vec<String>>()
        .join("\n")
}


use itertools::Itertools;
/// eijl2m (eijlist_to_matrix) takes a lexigraphically sorted list of transformations and returns the corresponding matrix representation
/// 
fn eijl2m<const K:usize>(eijlist: &Vec<(usize,usize)> ) -> [u32; K] {
    let mut res = [0; K];
    eijlist.iter().counts().into_iter().for_each(|(k,v)|
        { res[ (k.0 as f64).sqrt() as usize + k.1] = v as u32; }
    );
    res
}



#[derive(PartialEq, Eq, Clone, Debug)]
/// EijList is wrapper type for vector of (usize, usize) that allows it to be used as key in trie, automatically sorts it to use suffixes
struct EijList( Vec<(usize,usize)> );

impl TrieKey for EijList {
    fn encode_bytes(&self) -> Vec<u8> {
    let mut v = Vec::<u8>::with_capacity(self.0.len() * 2 * std::mem::size_of::<usize>());
    for u in self.0.iter().rev() { // Going in reverse here means we compare suffixes!
        v.extend_from_slice(&u.0.to_be_bytes());
        v.extend_from_slice(&u.1.to_be_bytes());
    }
    v
    }
}

/// apply_eij_prod_on_metavar applies the group action of a list of Eij onto a metavariable, by repeatedly calling apply_eij_on_metavar with short circuiting
///  
/// Note that the eij_prod should be passed smallest to largest lexigraphically as the function automatically applies .rev() to the iterator.
fn apply_eij_prod_on_metavar(
    m: &MetaVar,
    eijlist: &Vec<(usize, usize)>
) -> Option<(u32, MetaVar)> 
{
    eijlist.iter().rev().try_fold((1,m.clone()), |acc, eij| {
        apply_eij_on_metavar(&acc.1, *eij).map(|(x,m)| (x*acc.0,m))
    })
}

fn apply_eij_prod_on_input<T, F, const K: usize>(
    m: &MetaVar,
    lno: usize, // The line number of the input line
    eijlist: &Vec<(usize, usize)>,
    p_slp: &mut PartialSLP<T,K>,
    line_map: &mut LineMap<K>,
    i32_to_c: F,
)
where
    T: Clone + Display,
    F: Fn(i32) -> T + Copy,
{
    use Operation::*;
    use PartialSLPLine::*;
    use PartialSLPVar::*;
    if let Some((coeff, m_new)) = apply_eij_prod_on_metavar(m, eijlist) {
        // result is of form coeff * m_new, expressing in SLP is =m_new \n =Coeff*L(before)
        // partial_slp is constructed backwards meaning we add as follows: =Coeff*L(current slp len+1), =m_new
        //  the +1 is necessary as upon pushing, current slp len refers to index of line to be added, we want to refer to the line afterwards
        let p_slp_len = p_slp.len();
        p_slp.push(Compound((C(i32_to_c(coeff as i32)), Mult, LIP(p_slp_len + 1))));
        p_slp.push(Input(m_new));
        // add reference to the new entry point for Eij . m into the line mapper
        insert_lm::<K>(line_map, lno, eijl2m::<K>(&eijlist), p_slp_len);
    } else {
        // result is 0
        // know that by rules of SLP, a line cannot refer to itself => in p_slp will never have variable of LIP(0) hence this is reserved for 0
        insert_lm::<K>(line_map, lno, eijl2m::<K>(&eijlist), 0);
    }
}


// TODO:
// - actually use the lines_to_transform
// - actually use the tries hashmap to use cached results
// - change PartialSLP to allow references from the current slp_res (difficult)
// - implement partial_slp support for addition and multiplication
// - implement the conversion from partial_slp to slp

pub fn apply_eij_poly_on_program<T, F, const K: usize>(
    slp: &SLP<T>,
    eijs: &HashMap<Vec<(usize, usize)>, i32>,
    i32_to_c: F,
) -> Result<SLP<T>, String>
where
    T: Clone + Display,
    F: Fn(i32) -> T + Copy,
{
    use SLPLine::*;
    use SLPVar::*;
    use Operation::*;

    // create hashmap of line number and radix_trie of the combinations
    let mut tries: HashMap<usize, Trie<EijList, usize>> = HashMap::new();
    let mut slp_res:SLP<T> = vec![];

    'outer: for (prod,coeff) in eijs.iter().sorted_by(|(prod,coeff)| -prod.len()) { // go through each term in order of largest product to smallest
        let mut p_slp: PartialSLP<T,K> = vec![];
        let eijlist = EijList(prod.clone());

        let mut line_map: LineMap<K> = HashMap::new();
        let mut lines_to_transform: LinesToTransform<K> = HashMap::new();
        lines_to_transform.insert(slp.len()-1, eijl2m(&prod));

        // Using the trie:
        //  FOR THE VERY FIRST (LAST?) OF PROJECTING WITH THIS EIJLIST
        //    If we have already computed the entire thing, skip to adding that line number to the overall slp
        //  FOR EVERY OTHER LINE
        //    If the line is a Compound with references:
        //     - Find nearest ancestor for reference in trie
        //     - If ancestor key is eijlist, use the direct slp reference
        //     - If ancestor is suffix and is all but one (ie only on eij not in suffix) then can use basis_element code to compute from slp_res
        //     - If ancestor is suffix, can technically use from that line in slp_res (maybe extension) to skip further computation
        //     - If no ancestor then just have to compute it from scratch 
        for (lno, line) in slp.iter().enumerate().rev() {
            if !lines_to_transform.contains_key(&lno) { continue; }

            if lno == slp.len()-1 && let Some(trie) = tries.get(&lno) && let Some(&trans_lno) = trie.get(&eijlist) { // First line of computation, found reference in slp_res already
                slp_res.push( Compound(( C(i32_to_c(*coeff as i32)), Mult, L(trans_lno)       )) );
                slp_res.push( Compound(( L(slp_res.len()-2),         Plus, L(slp_res.len()-1) )) );
                continue 'outer; // already done for this product, move to next
            }

            if prod.is_empty() {
                
            }


            match line {
                Input(m) => apply_eij_prod_on_input(m, lno, prod, &mut p_slp, &mut line_map, i32_to_c), // TODO : ignoring the capability for trie
                Compound((C(_), Plus, C(_)) => (), // TODO : Continue from here
            }
        }



    }


    unimplemented!()
}
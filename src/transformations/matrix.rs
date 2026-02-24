//! matrix_transformations provides the singular operation of applying a linear combination of lists of sorted eijs to a slp.
//! This task is considered large enough to dedicate a file to ensure clarity.

#![allow(dead_code)]

use crate::{straight_line_program::*, transformations::apply_eij_on_metavar};

use radix_trie::Trie;
use std::collections::HashMap; 
use std::collections::HashSet;
use radix_trie::TrieKey;
use std::fmt::Display;

/// LinesToTransform holds dict of lines in (original / currently computing) SLP and which transforms need to be applied to them
type LinesToTransform<const K: usize> = HashMap<usize, HashSet<[u32;K]>>;
/// insert_l2t safely inserts elements into the LinesToTransform structure such that the transform is added to the set present
/// 
/// K: is the number of elements in the matrix (should be a square number K=k*k).
/// 
/// lm is the linemap.
/// ilno is the input slp line number.
/// transform is the transforming matrix.
fn insert_l2t<const K: usize>(l2t: &mut LinesToTransform<K>, ilno: usize, transform: [u32;K]) {
    if let Some(set) = l2t.get_mut(&ilno) {
        set.insert(transform);
    } else {
        l2t.insert(ilno, HashSet::from([transform]) );
    }
}
/// LineMap holds dict of lines and transformation applied from original SLP and location in new transformed partial SLP
type LineMap<const K: usize> = HashMap<usize, HashMap<[u32;K],usize>>;
/// insert_lm safely inserts elements into the LineMap structure such that the output line no. (olno) is inserted into the inner map
/// 
/// K: is the number of elements in the matrix (should be a square number K=k*k).
/// 
/// lm is the linemap.
/// ilno is the input slp line number.
/// transform is the transforming matrix.
/// olno is the output partial_slp line number to insert.
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
    LICP(usize), // LineInCurrentSLP is an index of an earlier line in the current slp the partial slp will be added to
    LTT((usize, [u32;K])), // LineToTranslate, [u32;K] represents the transformation
    LTTCP((usize, [u32;K])), // LineToTranslate_in_CurrentProgram, [u32;K] represents the transformation, usize the index
    C(T), // C(T) is the coefficient of type T
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
        LIP(lno) => format!("LIP{lno}"),
        LICP(lno) => format!("LICP{lno}"),
        LTT((lno, transform)) => format!("LTT({lno},T:{transform:?})"),
        LTTCP((lno, transform)) => format!("LTTCP({lno},T:{transform:?})"),
        C(c) => c.to_string(),
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
/// eijlist is the vector of transformations to be converted to matrix form
fn eijl2m<const K:usize>(eijlist: &Vec<(usize,usize)> ) -> [u32; K] {
    let mut res = [0; K];
    eijlist.iter().counts().into_iter().for_each(|(k,v)|
        { res[ (k.0 as f64).sqrt() as usize + k.1] = v as u32; }
    );
    res
}

/// eijm2l (eijmatrix_to_list) takes a k*k matrix and returns the corresponding  list of eij basis transformations
/// 
/// matrxi is the k*k matrix to convert to list form
fn eijm2l<const K:usize>(matrix: &[u32;K]) -> Vec<(usize,usize)> {
    let mut res = Vec::new();
    let k: usize = (K as f64).sqrt() as usize;
    matrix.iter().enumerate().for_each(|(i,j)|{
        [0..*j].iter().for_each(|_| res.push( (i/k, i%k) ));
    });
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

/// matrix_remove_last_n_transforms is a quick function to remove transformations from the end of the matrix
fn matrix_remove_last_n_transforms<const K: usize>(m: &mut [u32;K], n: u32) {
    for el in m.iter_mut().rev() {
        if *el > 0 {
            *el -= min(n, *el);
            n -= min(n, *el);
        }
        if n == 0 {
            break;
        }
    }
}

/// transform_lno_efficiently is some common functionality abstracted: it takes a line number (from a reference) that needs to be transformed
///  and finds the fastest way to get it transformed, including the trie to reduce the number of computations
/// 
/// l2t_input is the lines_to_transform where lines come from the original slp.
/// l2t_slp is the lines_to_transform where the lines come from the current slp it will be added to.
/// lno is the line number of the reference.
/// transform is the transform to apply.
/// tries is the map of the lines (from input slp) to the trie of transformations applied to it.
/// is_input_slp is boolean representing whether transforming line from the input slp or the current one
/// 
/// Returns the partial_slp variable that is the correct LineToTransform to transform the current line
fn transform_lno_efficiently<const K:usize>(
    l2t_input: &mut LinesToTransform<K>,
    l2t_slp: &mut LinesToTransform<K>,
    lno: usize,
    transform: [u32;K],
    tries: &HashMap<usize, Trie<EijList, usize>>,
    is_input_slp: bool,
) -> PartialSLPVar<K> {
    use PartialSLPVar::*;
    let eijlist = EijList( eijm2l::<K>(matrix) );
    if let Some(trie) = tries.get(&lno) && let Some(anc) = trie.get_ancestor( &eijlist ) {
        matrix_remove_last_n_transforms( &mut transform, anc.key().unwrap().len() );
        insert_l2t(l2t_slp, anc.value().unwrap(), transform.clone());

        LTTCP((anc.value().unwrap(), transform))
    } else {
        insert_l2t(l2t_input, lno, transform);

        if is_input_slp { LTT((lno, transform)) } else { LTTCP((lno, transform)) }
    }
}

/// apply_eij_on_addition computes the transformation of an addition line (=L2+0.5) under eij
/// This involves updating the partial structures used to build the final slp: p_slp and line_map
/// v1 the the first term in the line
/// v2 the the second term in the line
/// lno is the line number of this Input line in the original SLP, needed to give the key for line_map
/// p_slp is the partial slp structure being built to be converted into a slp, it is being built backwards
/// lines_to_trans is a dictionary of other lines that also need to be included into the partial slp
/// line_map is a dictionary mapping original line numbers and to their current index in the partial slp structure
fn apply_transform_on_addition<T, const K: usize>(
    v1: &SLPVar<T,K>,
    v2: &SLPVar<T,K>,
    transform: &[u32;K],
    lno: usize,
    p_slp: &mut PartialSLP<T,K>,
    l2t_input: &mut LinesToTransform<K>, // both should be l2t_slp for when transforming the current slp
    l2t_slp: &mut LinesToTransform<K>,
    curr_prog_to_trans: &mut LinesToTransform<K>,
    line_map: &mut LineMap<K>,
    tries: &HashMap<usize, Trie<EijList, usize>>,
    is_input_slp: bool,
) { // TODO : Write this logic on paper, ensures works
    use Operation::*;
    use PartialSLPLine::*;
    use PartialSLPVar::{LTT,LTTCP};
    use SLPVar::*;
    match (v1, v2) {
        (L(n1), L(n2)) => {
            let n1_var = transform_lno_efficiently(l2t_input, l2t_slp, *n1, transform.clone(), tries, is_input_slp);
            let n2_var = transform_lno_efficiently(l2t_input, l2t_slp, *n2, transform.clone(), tries, is_input_slp);
            let n1_var;
            if let Some((n1_lno,anc_len)) = n1_res {
                let t_clone = transform.clone();
                matrix_remove_last_n_transforms(&mut t_clone, anc_len as u32);
                n1_var = LTTCP((n1_lno, t_clone));
            } else if is_input_slp {
                n1_var = LTTCP((*n1, transform.clone()));
            } else {
                n1_var = LTT((*n1, transform.clone()));
            }

            p_slp.push(Compound(( n1_var, Plus, n2_var )));
            insert_lm(line_map, lno, transform.clone(), p_slp.len() - 1);
        }
        (L(n), C(_)) | (C(_), L(n)) => {
            let n_var = transform_lno_efficiently(l2t_input, l2t_slp, *n, transform.clone(), tries, is_input_slp);
            p_slp.push(Compound(( n_var, Plus, LIP(0) )));

            insert_lm(line_map, lno, transform.clone(), p_slp.len() - 1);
        }
        (C(_), C(_)) => insert_lm(line_map, lno, transform.clone(), 0),
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
        lines_to_transform.insert(slp.len()-1, eijl2m::<K>(&prod));

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
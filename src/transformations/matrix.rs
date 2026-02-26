//! matrix_transformations provides the singular operation of applying a linear combination of lists of sorted eijs to a slp.
//! This task is considered large enough to dedicate a file to ensure clarity.

#![allow(dead_code)]

use crate::{straight_line_program::*, transformations::apply_eij_on_metavar};

use radix_trie::Trie;
use radix_trie::TrieCommon;
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

/// TrieMap maps an input line to a trie storing all transformations performed on that line alongside their index in the current program
type TrieMap = HashMap<usize, Trie<EijList, usize>>;

/// insert_tm inserts a value into the map of input lines to tries
fn insert_tm<const K: usize>(tm: &mut TrieMap, ilno: usize, transform: &[u32;K], olno: usize) {
    if let Some(trie) = tm.get_mut(&ilno) {
        if let Some(_) = trie.get( &EijList(eijm2l(transform)) ) { println!("Trie storage is trying to overwrite!!!");}
        trie.insert( EijList(eijm2l(transform)), olno);
    } else {
        let mut trie = Trie::new();
        trie.insert( EijList(eijm2l(transform)), olno);

        tm.insert(ilno, trie);
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
    let k = (K as f64).sqrt() as usize;
    eijlist.iter().counts().into_iter().for_each(|(p,v)|
        { res[ p.0*k + p.1] = v as u32; }
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
        (0..*j).for_each(|_| res.push( (i/k, i%k) ));
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
) -> Option<(u64, MetaVar)> 
{
    eijlist.iter().rev().try_fold((1 as u64,m.clone()), |acc, eij| {
        apply_eij_on_metavar(&acc.1, *eij).map(|(x,m)| ( x as u64 *acc.0 ,m))
    })
}

fn apply_eij_prod_on_input<T, F, const K: usize>(
    m: &MetaVar,
    lno: usize, // The line number of the input line
    eijlist: &Vec<(usize, usize)>,
    p_slp: &mut PartialSLP<T,K>,
    line_map: &mut LineMap<K>,
    u64_to_c: F,
)
where
    T: Clone + Display,
    F: Fn(u64) -> T + Copy,
{
    use Operation::*;
    use PartialSLPLine::*;
    use PartialSLPVar::*;
    if let Some((coeff, m_new)) = apply_eij_prod_on_metavar(m, eijlist) {
        // result is of form coeff * m_new, expressing in SLP is =m_new \n =Coeff*L(before)
        // partial_slp is constructed backwards meaning we add as follows: =Coeff*L(current slp len+1), =m_new
        //  the +1 is necessary as upon pushing, current slp len refers to index of line to be added, we want to refer to the line afterwards
        let p_slp_len = p_slp.len();
        p_slp.push(Compound((C(u64_to_c(coeff)), Mult, LIP(p_slp_len + 1))));
        p_slp.push(Input(m_new));
        // add reference to the new entry point for Eij . m into the line mapper
        insert_lm(line_map, lno, eijl2m::<K>(&eijlist), p_slp_len);
    } else {
        // result is 0
        // know that by rules of SLP, a line cannot refer to itself => in p_slp will never have variable of LIP(0) hence this is reserved for 0
        insert_lm(line_map, lno, eijl2m::<K>(&eijlist), 0);
    }
}

/// matrix_remove_last_n_transforms is a quick function to remove transformations from the end of the matrix
fn matrix_remove_last_n_transforms<const K: usize>(m: &mut [u32;K], mut n: u32) {
    for el in m.iter_mut().rev() {
        if *el > 0 {
            let min = n.min(*el);
            *el -= min;
            n -= min;
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
/// l2t_slp is the option dictionary of other lines in the input program that also need to be included into the partial slp with their correspnding transform. If none represents that l2t_input is l2t_slp
/// lno is the line number of the reference.
/// transform is the transform to apply.
/// tries is the map of the lines (from input slp) to the trie of transformations applied to it.
/// 
/// Returns the partial_slp variable that is the correct LineToTransform to transform the current line
fn transform_lno_efficiently<T: Clone, const K:usize>(
    l2t_input: &mut LinesToTransform<K>, // if processing lines in the current slp, use l2t_slp here, leave l2t_slp as None
    l2t_slp: &mut Option<&mut LinesToTransform<K>>,
    lno: usize,
    transform: [u32;K],
    tries: &HashMap<usize, Trie<EijList, usize>>,
) -> PartialSLPVar<T,K> {
    use PartialSLPVar::*;
    let eijlist = EijList( eijm2l::<K>(&transform) );
    if let Some(trie) = tries.get(&lno) && let Some(clno) = trie.get(&eijlist) { // have literally computed this line + transform before!!
        LICP(*clno)
    } else if let Some(trie) = tries.get(&lno) && let Some(anc) = trie.get_ancestor( &eijlist ) { // have computed this line with subtransform before
        let mut t_clone = transform.clone();
        matrix_remove_last_n_transforms( &mut t_clone, anc.key().unwrap().0.len() as u32 );
        if let Some(l2t_slp) = l2t_slp.as_mut() {
            insert_l2t(l2t_slp, *anc.value().unwrap(), t_clone.clone());
        } else {
            insert_l2t(l2t_input, *anc.value().unwrap(), t_clone.clone());
        }

        LTTCP((*anc.value().unwrap(), t_clone))
    } else {
        insert_l2t(l2t_input, lno, transform);

        if l2t_slp.is_some() { LTT((lno, transform)) } else { LTTCP((lno, transform)) }
    }
}

// TODO : Write this logic on paper, ensures works
/// apply_eij_prod_on_addition computes the transformation of an addition line (=L2+0.5) under a matrix representing sorted list of eij transformations
/// This involves updating the partial structures used to build the final slp: p_slp and line_map.
/// 
/// v1 the the first term in the line.
/// v2 the the second term in the line.
/// lno is the line number of this Input line in the original SLP, needed to give the key for line_map.
/// p_slp is the partial slp structure being built to be converted into a slp, it is being built backwards.
/// l2t_input is a dictionary of other lines in the input program that also need to be included into the partial slp with their correspnding transform.
/// l2t_slp is the option dictionary of other lines in the input program that also need to be included into the partial slp with their correspnding transform. If none represents that l2t_input is l2t_slp
/// line_map is a dictionary mapping original line numbers to their current index in the partial slp structure.
/// tries is a map of input line numbers to a trie of transforms and their output line in the current slp.
fn apply_eij_prod_on_addition<T: Clone, const K: usize>(
    v1: &SLPVar<T>,
    v2: &SLPVar<T>,
    transform: &[u32;K],
    lno: usize,
    p_slp: &mut PartialSLP<T,K>,
    l2t_input: &mut LinesToTransform<K>,
    l2t_slp: &mut Option<&mut LinesToTransform<K>>,
    line_map: &mut LineMap<K>,
    tries: &HashMap<usize, Trie<EijList, usize>>,
) { 
    use Operation::*;
    use PartialSLPLine::*;
    use PartialSLPVar::{LIP,LTT,LTTCP};
    use SLPVar::*;
    match (v1, v2) {
        (L(n1), L(n2)) => {
            let n1_var = transform_lno_efficiently(l2t_input, l2t_slp, *n1, transform.clone(), tries);
            let n2_var = transform_lno_efficiently(l2t_input, l2t_slp, *n2, transform.clone(), tries);

            p_slp.push(Compound(( n1_var, Plus, n2_var )));
            insert_lm(line_map, lno, transform.clone(), p_slp.len() - 1);
        }
        (L(n), C(_)) | (C(_), L(n)) => {
            let n_var = transform_lno_efficiently(l2t_input, l2t_slp, *n, transform.clone(), tries);
            p_slp.push(Compound(( n_var, Plus, LIP(0) )));

            insert_lm(line_map, lno, transform.clone(), p_slp.len() - 1);
        }
        (C(_), C(_)) => insert_lm(line_map, lno, transform.clone(), 0),
    }
}

/// get_all_sub_transforms takes a matrix and returns all others with entries smaller than the original
fn get_all_sub_transforms<const K:usize>(m: [u32;K]) -> Vec<[u32;K]> {
    let mut q: Vec<[u32;K]> = vec![m];
    for (i,el) in m.iter().enumerate() {
        if *el == 0 { continue; }
        q = q.iter().map(|v|
            (0..=*el).map(|j| {
                let mut clone = v.clone();
                clone[i] = j;
                clone
            })
        ).flatten().collect();
    }
    q
}

/// apply_eij_prod_on_mult computes the transformation of a multiplication line (=L2*0.5) under a matrix representing sorted list of eij transformations
/// This involves updating the partial structures used to build the final slp: p_slp and line_map.
/// v1 the the first term in the line.
/// v2 the the second term in the line.
/// lno is the line number of this Input line in the original SLP, needed to give the key for line_map.
/// p_slp is the partial slp structure being built to be converted into a slp, it is being built backwards.
/// l2t_input is a dictionary of other lines in the input program that also need to be included into the partial slp with their correspnding transform.
/// l2t_slp is a dictionary of other lines in the input program that also need to be included into the partial slp with their correspnding transform.
/// line_map is a dictionary mapping original line numbers to their current index in the partial slp structure.
/// tries is a map of input line numbers to a trie of transforms and their output line in the current slp.
fn apply_eij_prod_on_mult<T: Clone, const K: usize>(
    v1: &SLPVar<T>,
    v2: &SLPVar<T>,
    transform: &[u32;K],
    lno: usize,
    p_slp: &mut PartialSLP<T,K>,
    l2t_input: &mut LinesToTransform<K>, // both should be l2t_slp for when transforming the current slp
    l2t_slp: &mut Option<&mut LinesToTransform<K>>,
    line_map: &mut LineMap<K>,
    tries: &HashMap<usize, Trie<EijList, usize>>,
) {
    use Operation::*;
    use PartialSLPLine::*;
    use PartialSLPVar::C as PC;
    use PartialSLPVar::{LIP,LTT};
    use SLPVar::*;
    match (v1, v2) {
        (L(n1), L(n2)) => {
            let all_subvecs = get_all_sub_transforms(transform.clone());
            let all_subvecs_len = all_subvecs.len();
            for (i, subvec) in all_subvecs.into_iter().enumerate() {
                let mut subvec_minus = [0;K];
                subvec_minus.iter_mut().enumerate().for_each(|(i,el)| *el = transform[i]-subvec[i]);

                let n1_var = transform_lno_efficiently(l2t_input, l2t_slp, *n1, subvec, tries);
                let n2_var = transform_lno_efficiently(l2t_input, l2t_slp, *n2, subvec_minus, tries);
                if i != all_subvecs_len-1 { // if not the final iteration
                    p_slp.push(Compound(( LIP(p_slp.len()+1), Plus, LIP(p_slp.len()+2) )));
                }
                p_slp.push(Compound(( n1_var, Mult, n2_var )));

            }

            insert_lm(line_map, lno, transform.clone(), p_slp.len()-1);
        }
        (L(n), C(c)) | (C(c), L(n)) => {
            let n_var = transform_lno_efficiently(l2t_input, l2t_slp, *n, transform.clone(), tries);
            p_slp.push(Compound(( PC(c.clone()), Mult, n_var )));

            insert_lm(line_map, lno, transform.clone(), p_slp.len()-1);
        }
        (C(_), C(_)) => insert_lm(line_map, lno, transform.clone(), 0),
    }
}

/// apply_eij_prod_on_line computes the transformation of a general line (=L2*0.5) under a matrix representing sorted list of eij transformations
/// This involves updating the partial structures used to build the final slp: p_slp and line_map.
/// 
/// v1 the the first term in the line.
/// v2 the the second term in the line.
/// lno is the line number of this Input line in the original SLP, needed to give the key for line_map.
/// p_slp is the partial slp structure being built to be converted into a slp, it is being built backwards.
/// l2t_input is a dictionary of other lines in the input program that also need to be included into the partial slp with their correspnding transform.
/// l2t_slp is a dictionary of other lines in the input program that also need to be included into the partial slp with their correspnding transform.
/// line_map is a dictionary mapping original line numbers to their current index in the partial slp structure.
/// tries is a map of input line numbers to a trie of transforms and their output line in the current slp.
/// is_input_slp is a boolean value representing whether the line comes from the input slp or the currently constructed one.
/// i32_to_c is a function that converts a i32 value to a value in your chosen type, it should respect that i32_to_c(0) = 0.
/// 
/// Note that the transform used should not be 0_{n*n} as this simpler case is not handled by the functions this function calls (this transform should instead be separately handled).
fn apply_eij_prod_on_line<T, F, const K: usize>(
    line: &SLPLine<T>,
    transform: &[u32;K],
    lno: usize,
    p_slp: &mut PartialSLP<T,K>,
    l2t_input: &mut LinesToTransform<K>, // both should be l2t_slp for when transforming the current slp
    l2t_slp: &mut Option<&mut LinesToTransform<K>>,
    line_map: &mut LineMap<K>,
    tries: &HashMap<usize, Trie<EijList, usize>>,
    u64_to_c: F,
)
where T: Clone + Display, F: Fn(u64)->T + Copy,
{
    use SLPLine::*;
    use Operation::*;
    match line {
        Input(m) => apply_eij_prod_on_input(m, lno, &eijm2l(transform), p_slp, line_map, u64_to_c),
        Compound((v1,Plus,v2)) => apply_eij_prod_on_addition(v1, v2, transform, lno, p_slp, l2t_input, l2t_slp, line_map, tries),
        Compound((v1, Mult, v2)) => apply_eij_prod_on_mult(v1, v2, transform, lno, p_slp, l2t_input, l2t_slp, line_map, tries),
    }
}

/// include_line_from_slp includes a line from the SLP as is into the partial slp
/// This means copying the line across and ensuring other lines referenced will be copied across also
/// line the the line being computed
/// lno is the line number of this Input line in the original SLP, needed to give the key for line_map
/// p_slp is the partial slp structure being built to be converted into a slp, it is being built backwards
/// lines_to_trans is a dictionary of other lines that also need to be included into the partial slp
/// line_map is a dictionary mapping original line numbers and to their current index in the partial slp structure
fn include_line_from_slp<T: Clone, const K: usize>(
    line: &SLPLine<T>,
    transform: &[u32;K],
    lno: usize,
    p_slp: &mut PartialSLP<T,K>,
    l2t_input: &mut LinesToTransform<K>, // both should be l2t_slp for when transforming the current slp
    l2t_slp: &mut Option<&mut LinesToTransform<K>>,
    line_map: &mut LineMap<K>,
    tries: &HashMap<usize, Trie<EijList, usize>>,
) {
    use PartialSLPLine::Compound as PCompound;
    use PartialSLPLine::Input as PInput;
    use PartialSLPVar::C as PC;
    use PartialSLPVar::LTT;
    use SLPLine::*;
    use SLPVar::*;
    match line {
        Compound((L(n1), op, L(n2))) => {
            let n1_var = transform_lno_efficiently(l2t_input, l2t_slp, *n1, transform.clone(), tries);
            let n2_var = transform_lno_efficiently(l2t_input, l2t_slp, *n2, transform.clone(), tries);
            p_slp.push(PCompound(( n1_var, *op, n2_var )));
        }
        Compound((L(n), op, C(c))) | Compound((C(c), op, L(n))) => {
            let n_var = transform_lno_efficiently(l2t_input, l2t_slp, *n, transform.clone(), tries);
            p_slp.push(PCompound((PC(c.clone()), *op, n_var )));
        }
        Compound((C(c1), op, C(c2))) => {
            p_slp.push(PCompound((PC(c1.clone()), *op, PC(c2.clone()))))
        }
        Input(m) => p_slp.push(PInput(m.clone())),
    }
    insert_lm(line_map, lno, transform.clone(), p_slp.len()-1);
}

fn apply_eij_prod_on_slp<T, F, const K: usize>(
    slp: &SLP<T>,
    curr_slp: &SLP<T>,
    transform: &[u32;K],
    tries: &HashMap<usize, Trie<EijList, usize>>,
    u64_to_c: F
) -> Result< (PartialSLP<T,K>, LinesToTransform<K>, LinesToTransform<K>, LineMap<K>, LineMap<K>), String>
where
    T: Clone + Display,
    F: Fn(u64) -> T + Copy,
{

    // Steps to convert SLP:
    //  0. Validity of inputs
    //  1. Initialisation of structures
    //  2. Calculation of each line in input slp
    //  3. Calculation of each line needed in current slp
    //  4. Return the structures holding all components of the current slp together

    // 0. Validity of inputs
    if let Some(i) = are_all_line_refs_valid(slp) {
        return Err(format!(
            "SLP contains line reference to a line not strictly before itself on line {i}:{}",
            stringify_slpline(&slp[i])
        ));
    }
    if !are_metavars_all_same_len(slp) {
        return Err("SLP contains metavariables of different lengths".into());
    }
    if !do_metavars_refer_to_homogeneous_poly(slp) {
        return Err("SLP has metavariables referring to non-homogeneous polynomial".into());
    }
    let k = (K as f64).sqrt() as usize;
    if k != get_metavar_len(slp).unwrap() {
        return Err(format!("Transform is {k}*{k} matrix but metavariables use {} variables instead", get_metavar_len(slp).unwrap()) );
    }

    // 1. Initialisation of structures
    let mut p_slp: PartialSLP<T,K> = vec![];
    let mut l2t_input: LinesToTransform<K> = HashMap::new();
    let mut l2t_slp: LinesToTransform<K> = HashMap::new();
    let mut input_line_map: LineMap<K> = HashMap::new();
    let mut curr_line_map: LineMap<K> = HashMap::new();

    insert_l2t(&mut l2t_input, slp.len()-1, transform.clone() );

    //  2. Calculation of each line in input slp
    for (lno, line) in slp.iter().enumerate().rev() {
        for transform in l2t_input.get(&lno).cloned().into_iter().flatten() { // edit here to change the order that transforms are applied in
            if transform.iter().all(|i| *i==0) {
                include_line_from_slp( line, &transform, lno, &mut p_slp, &mut l2t_input, &mut Some(&mut l2t_slp), &mut input_line_map, tries);
            } else {
                apply_eij_prod_on_line(line, &transform, lno, &mut p_slp, &mut l2t_input, &mut Some(&mut l2t_slp), &mut input_line_map, tries, u64_to_c);
            }
        }
    }

    //3. Calculation of each line needed in current slp
    for (lno, line) in curr_slp.iter().enumerate().rev() {
        for transform in l2t_slp.get(&lno).cloned().into_iter().flatten() { // edit here to change the order that transforms are applied in
            if transform.iter().all(|i| *i==0) {
                include_line_from_slp( line, &transform, lno, &mut p_slp, &mut l2t_slp, &mut None, &mut curr_line_map, tries);
            } else {
                apply_eij_prod_on_line(line, &transform, lno, &mut p_slp, &mut l2t_slp, &mut None, &mut curr_line_map, tries, u64_to_c);
            }
        }
    }

    Ok((p_slp, l2t_input, l2t_slp, input_line_map, curr_line_map))
}

fn add_scaled_p_slp_onto_curr_slp<T, const K: usize>(
    mut curr_slp: SLP<T>,
    mut p_slp: PartialSLP<T,K>,
    scale: T,
    // l2t_input: LinesToTransform<K>,
    // l2t_slp: LinesToTransform<K>,
    input_line_map: &mut LineMap<K>,
    curr_line_map: &LineMap<K>,
    tries: &mut TrieMap,
    zero_coeff: T,
) -> Result<SLP<T>, String>
where
    T: Clone + Display,
{
    use PartialSLPLine::*;
    use PartialSLPVar::*;
    use Operation::*;

    // 3. Inversion and conversion
    p_slp.push(Compound(( C(zero_coeff.clone()), Plus, C(zero_coeff.clone()) ))); // add zero to the beginning (once reversed)
    p_slp.reverse();

    // convert all LTT to LIP, LTTCP to LIP as at this point every line should have a transformation in the p_slp
    for line in p_slp.iter_mut() {
        match line {
            Compound((v1, _, v2)) => {
                match v1 {
                    LTT(  (p1, m)) => *v1 = LIP(*input_line_map.get(p1).ok_or("Couldn't find transformation")?.get(m).ok_or("Couldn't find transformation")?),
                    LTTCP((p1, m)) => *v1 = LIP( *curr_line_map.get(p1).ok_or("Couldn't find transformation")?.get(m).ok_or("Couldn't find transformation")?),
                    _ => (), 
                }
                match v2 {
                    LTT(  (p2, m)) => *v1 = LIP(*input_line_map.get(p2).ok_or("Couldn't find transformation")?.get(m).ok_or("Couldn't find transformation")?),
                    LTTCP((p2, m)) => *v1 = LIP( *curr_line_map.get(p2).ok_or("Couldn't find transformation")?.get(m).ok_or("Couldn't find transformation")?),
                    _ => (), 
                }
            }
            Input(_) => (),
        }
    }

    // all LIP in the p_slp need to be reversed to point back to where they should
    // all LIP need to be increased by one as well as referred to location in SLP without added 0 element
    //  1 -> p_len-2 and p_len-1 -> 0, note that 0 shouldn't be reversed as it points to an element that is now at the beginning
    //  so use equation f(x) = p_slp.len() - x - 1
    let p_len = p_slp.len();
    for line in p_slp.iter_mut() {
        match line {
            Compound((LIP(0), _, LIP(0))) => (), // handle cases where we don't want to map 0 first
            Compound((LIP(0), _, LIP(n))) | Compound((LIP(n), _, LIP(0))) => *n = p_len - *n - 1,
            Compound((LIP(0), _, _)) | Compound((_, _, LIP(0))) => (),
            Compound((LIP(n1), _, LIP(n2))) => {
                *n1 = p_len - *n1 - 1;
                *n2 = p_len - *n2 - 1;
            }
            Compound((LIP(n), _, _)) | Compound((_, _, LIP(n))) => *n = p_len - *n - 1,
            _ => (),
        }
    }
    // also will map the references in the line maps enabling addition of these structure to tries
    input_line_map.iter_mut().for_each(|(_k,v)| 
        v.values_mut().filter(|lno| **lno!=0).for_each(|lno|
            *lno = curr_slp.len() + p_len - *lno - 1
        )
    );

    for (ilno, tms) in input_line_map {
        for (tm, olno) in tms {
            insert_tm(tries, *ilno, &tm, *olno);
        }
    }

    // convert PartialSLP to SLP, should be case that all variables are LIP, LICP or C
    use SLPLine::Compound as SLPCompound;
    use SLPLine::Input as SLPInput;
    use SLPVar::C as SLPC;
    use SLPVar::L;
    let mut new_slp: SLP<T> = Vec::new();
    let pv_to_sv = |v: PartialSLPVar<T, K>| match v {
        LIP(n) => Ok( L(n + curr_slp.len()) ),
        LICP(n) => Ok(L(n)),
        C(c) => Ok(SLPC(c.clone())),
        LTT(_) => Err("Unconverted LTT found"),
        LTTCP(_) => Err("Unconverted LTTCP found"),
    };
    for line in p_slp {
        new_slp.push(match line {
            Compound((v1, op, v2)) => SLPCompound((pv_to_sv(v1)?, op, pv_to_sv(v2)?)),
            Input(m) => SLPInput(m.clone()),
        });
    }

    let curr_len = curr_slp.len();
    curr_slp.append( &mut new_slp );
    curr_slp.push( SLPCompound(( SLPC(scale), Mult, L(curr_slp.len()-1) )));
    curr_slp.push( SLPCompound(( L(curr_len-1), Plus, L(curr_slp.len()-1) )));

    Ok(curr_slp)
}

// TODO: TEST THE APPLY_EIJ_PROD_ON_... FUNCTIONS!!
// TODO:
// - actually use the lines_to_transform
// - actually use the tries hashmap to use cached results
// - change PartialSLP to allow references from the current slp_res (difficult)
// - implement partial_slp support for addition and multiplication
// - implement the conversion from partial_slp to slp


// pub fn apply_eij_poly_on_program<T, F, const K: usize>(
//     slp: &SLP<T>,
//     eijs: &HashMap<Vec<(usize, usize)>, i32>,
//     i32_to_c: F,
// ) -> Result<SLP<T>, String>
// where
//     T: Clone + Display,
//     F: Fn(i32) -> T + Copy,
// {
//     use SLPLine::*;
//     use SLPVar::*;
//     use Operation::*;

//     // create hashmap of line number and radix_trie of the combinations
//     let mut tries: HashMap<usize, Trie<EijList, usize>> = HashMap::new();
//     let mut slp_res:SLP<T> = vec![];

//     'outer: for (prod,coeff) in eijs.iter().sorted_by(|(prod,coeff)| -prod.len()) { // go through each term in order of largest product to smallest
//         let mut p_slp: PartialSLP<T,K> = vec![];
//         let eijlist = EijList(prod.clone());

//         let mut line_map: LineMap<K> = HashMap::new();
//         let mut lines_to_transform: LinesToTransform<K> = HashMap::new();
//         lines_to_transform.insert(slp.len()-1, eijl2m::<K>(&prod));

//         // Using the trie:
//         //  FOR THE VERY FIRST (LAST?) OF PROJECTING WITH THIS EIJLIST
//         //    If we have already computed the entire thing, skip to adding that line number to the overall slp
//         //  FOR EVERY OTHER LINE
//         //    If the line is a Compound with references:
//         //     - Find nearest ancestor for reference in trie
//         //     - If ancestor key is eijlist, use the direct slp reference
//         //     - If ancestor is suffix and is all but one (ie only on eij not in suffix) then can use basis_element code to compute from slp_res
//         //     - If ancestor is suffix, can technically use from that line in slp_res (maybe extension) to skip further computation
//         //     - If no ancestor then just have to compute it from scratch 
//         for (lno, line) in slp.iter().enumerate().rev() {
//             if !lines_to_transform.contains_key(&lno) { continue; }

//             if lno == slp.len()-1 && let Some(trie) = tries.get(&lno) && let Some(&trans_lno) = trie.get(&eijlist) { // First line of computation, found reference in slp_res already
//                 slp_res.push( Compound(( C(i32_to_c(*coeff as i32)), Mult, L(trans_lno)       )) );
//                 slp_res.push( Compound(( L(slp_res.len()-2),         Plus, L(slp_res.len()-1) )) );
//                 continue 'outer; // already done for this product, move to next
//             }

//             if prod.is_empty() {
                
//             }


//             match line {
//                 Input(m) => apply_eij_prod_on_input(m, lno, prod, &mut p_slp, &mut line_map, i32_to_c), // TODO : ignoring the capability for trie
//                 Compound((C(_), Plus, C(_)) => (), // TODO : Continue from here
//             }
//         }



//     }


//     unimplemented!()
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evaluation;
    use crate::parsing;
    use num_rational::Rational64;

    use SLPLine::*;
    use SLPVar::*;

    #[test]
    fn eij_trie_test() {
        let mut trie: Trie<EijList, usize> = Trie::new();
        let vec1 = EijList( vec![(0,0), (0,0), (1,1), (1,2)] );
        let vec2 = EijList( vec![(1,1), (1,2)] );
        let vec3 = EijList( vec![(1,2)] );
        trie.insert(vec2.clone(), 2);
        trie.insert(vec3.clone(), 3);

        println!("nearest ancestor of {vec1:?} is {:?}",trie.get_ancestor(&vec1).unwrap().key());
        println!("nearest ancestor of {vec2:?} is {:?}",trie.get_ancestor(&vec2).unwrap().key());
    }

    #[test]
    fn eij_prod_input_test() {
        let metavar =  vec![50,50,50];
        let transform = [0,0,2,2,0,0,0,0,0];
        println!("eijlist: {:?}", eijm2l(&transform));
        println!("eij matrix {:?}", eijl2m::<9>( &eijm2l(&transform) ));
        let mut p_slp = vec![];
        let mut line_map: LineMap<9> = HashMap::new();
        let u64_to_c = |i| i as u64;
        apply_eij_prod_on_input::<u64, _, 9>(&metavar, 1, &eijm2l(&transform), &mut p_slp, &mut line_map, u64_to_c);

        println!("p_slp: {}", stringify_partialslp(&p_slp));
        println!("line map: {:?}", line_map);
    }

    #[test]
    fn eij_prod_addition_test() {
        let var1 = L(101);
        let var2 = C(104);
        let transform = [0,0,2,2,0,0,0,0,0];
        println!("eijlist: {:?}", eijm2l(&transform));
        println!("eij matrix {:?}", eijl2m::<9>( &eijm2l(&transform) ));

        let mut p_slp = vec![];

        let mut l2t_input: LinesToTransform<9> = HashMap::new();
        let mut l2t_slp: LinesToTransform<9> = HashMap::new();
        let mut line_map: LineMap<9> = HashMap::new();
        let mut trie_map: TrieMap = HashMap::new();
        insert_tm::<9>(&mut trie_map, 101, &[0,0,0,2,0,0,0,0,0], 20);

        apply_eij_prod_on_addition::<u64,9>(&var1, &var2, &transform, 10, &mut p_slp, &mut l2t_input, &mut Some(&mut l2t_slp), &mut line_map, &trie_map);

        println!("p_slp: {}", stringify_partialslp(&p_slp));
        println!("l2t_input: {l2t_input:?}");
        println!("l2t_slp: {l2t_slp:?}");
        println!("line map: {:?}", line_map);
    }

    #[test]
    fn eij_prod_mult_test() {
        let var1 = L(101);
        let var2 = L(104);
        let transform = [0,0,2,2,0,0,0,0,0];
        println!("eijlist: {:?}", eijm2l(&transform));
        println!("eij matrix {:?}", eijl2m::<9>( &eijm2l(&transform) ));

        let mut p_slp = vec![];

        let mut l2t_input: LinesToTransform<9> = HashMap::new();
        let mut l2t_slp: LinesToTransform<9> = HashMap::new();
        let mut line_map: LineMap<9> = HashMap::new();
        let mut trie_map: TrieMap = HashMap::new();
        insert_tm::<9>(&mut trie_map, 101, &[0,0,0,2,0,0,0,0,0], 20);

        apply_eij_prod_on_mult::<u64,9>(&var1, &var2, &transform, 10, &mut p_slp, &mut l2t_input, &mut Some(&mut l2t_slp), &mut line_map, &trie_map);

        println!("p_slp: {}", stringify_partialslp(&p_slp));
        println!("l2t_input: {l2t_input:?}");
        println!("l2t_slp: {l2t_slp:?}");
        println!("line map: {:?}", line_map);
    }



}
//! matrix_transformations provides the singular operation of applying a linear combination of lists of sorted eijs to a slp.
//! This task is considered large enough to dedicate a file to ensure clarity.

use crate::straight_line_program::*;
use super::basis_element::PartialSLP;

use radix_trie::Trie;
use std::collections::{HashMap, VecDeque};
use radix_trie::TrieKey;
use std::fmt::Display;


#[derive(PartialEq, Eq, Clone, Debug)]
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

pub fn apply_eij_poly_on_program<T, F>(
    slp: &SLP<T>,
    eijs: &HashMap<Vec<(usize, usize)>, i64>,
    i32_to_c: F,
) -> Result<SLP<T>, String>
where
    T: Clone + Display,
    F: Fn(i32) -> T + Copy,
{
    // create hashmap of line number and radix_trie of the combinations
    let mut tries:HashMap<usize, Trie<EijList, usize>> = HashMap::new();
    let mut slp_res:SLP<T> = vec![];

    for (prod,coeff) in eijs {
        let mut transform_q: VecDeque<(usize, EijList)> = VecDeque::new();
        transform_q.push_back(( slp.len()-1, EijList(prod.clone()) ));
        let mut p_slp: PartialSLP<T> = vec![];

        // Using the trie:
        //  FOR THE VERY FIRST (LAST?) OF PROJECTING WITH THIS EIJLIST
        //    If we have already computed the entire thing, skip to adding that line number to the overall slp
        //  FOR EVERY OTHER LINE
        //    If the line is a Compound with references:
        //     - Find nearest ancestor for reference in trie
        //     - If ancestor key is eijlist, use the direct slp reference
        //     - If ancestor is suffix, 
        while let Some((lno, eijlist)) = transform_q.pop_back() {
        }



    }


    unimplemented!()
}
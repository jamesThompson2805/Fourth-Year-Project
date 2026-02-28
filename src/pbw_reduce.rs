
use std::cmp::Ordering;
use std::collections::HashMap;
use itertools::Itertools;
use std::collections::BTreeSet;

/// SwapsItem allows the sorting of lists of swaps by their length, allowing the usage of the binary heap structure
///  items of same length are sorted lexigraphically, meaning items furthest to being sorted are sorted first and hence we won't repeat
#[derive(Eq,PartialEq,Clone,Hash,Debug)]
struct SwapsItem(Vec<(usize,usize)>);
impl Ord for SwapsItem {
    fn cmp(&self, other: &Self) -> Ordering {
        if other.0.len() == self.0.len() { // if lengths are same, reports it lexigraphically
            self.0.cmp(&other.0)
        } else {
            self.0.len().cmp(&other.0.len())
        }
    }
}
impl PartialOrd for SwapsItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// NOTE: This list of lists tells us how to accomplish partial computation as well, any suffix in this list implies 
//        we can compute that suffix in the slp first and appropriately reference the correct line, skipping calculation already done.

/// pbw_reduce takes an unsorted array of basis elements E_ij and sorts it, adding other sorted terms such that the transformations stay equal
/// swaps is the array of basis elements
/// returns hashmap representing the terms in the sum of an equation equal to the original swaps, all other terms are of smaller length 
pub fn pbw_reduce(swaps: &[(usize,usize)]) -> HashMap<Vec<(usize,usize)>, i32> {
    let mut bin_heap: BTreeSet<SwapsItem> =  BTreeSet::new();
    let mut signs: HashMap<Vec<(usize,usize)>, i32> =  HashMap::new();
    let mut result: HashMap<Vec<(usize,usize)>, i32> =  HashMap::new();
    bin_heap.insert( SwapsItem(Vec::from(swaps)) );
    signs.insert(Vec::from(swaps),1);

    while let Some(SwapsItem(mut swaps)) = bin_heap.pop_last() {
        let coeff = signs.remove(&swaps).unwrap();
        if coeff == 0 { // if coeff is 0 then this computation should be skipped
            continue;
        }
        // find first adjacent pair not lexigraphically in order by index
        let i = swaps.windows(2).find_position(|w| w[0] > w[1] ).map(|p| p.0);
        if let Some(i) = i { // found index needing swapping
            // ..,(i,j),(k,l),.. -> ..,(k,l),(i,j),.. + ..,[(i,j),(k,l)],..
            if swaps[i].1 == swaps[i+1].0 {
                let mut term1 =  swaps.clone();
                term1[i] = ( swaps[i].0, swaps[i+1].1 );
                term1.remove(i+1);
                if let Some(s) = signs.get_mut(&term1) {
                    *s += coeff;
                } else { signs.insert(term1.clone(), coeff); }
                bin_heap.insert(SwapsItem(term1));
            }
            if swaps[i].0 == swaps[i+1].1 {
                let mut term2 =  swaps.clone();
                term2[i] = ( swaps[i+1].0, swaps[i].1 );
                term2.remove(i+1);
                if let Some(s) = signs.get_mut(&term2) {
                    *s -= coeff;
                } else { signs.insert(term2.clone(), -coeff); }
                bin_heap.insert(SwapsItem(term2));
            }

            swaps.swap(i,i+1);
            if let Some(s) = signs.get_mut(&swaps) {
                *s += coeff;
            } else { signs.insert(swaps.clone(), coeff); }
            bin_heap.insert(SwapsItem(swaps));

        } else { // no index needing swapping: sorted this element
            result.insert(swaps, coeff);
        }
    }

    result
}
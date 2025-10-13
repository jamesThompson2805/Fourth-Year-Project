use crate::straight_line_program::MetaVar;
use crate::straight_line_program::LineVar;
use crate::straight_line_program::Operation::{Plus,Mult};
use crate::straight_line_program::SLPLine;
use crate::straight_line_program::SLPVar::{L,F};
use crate::straight_line_program::SLP;

use nom::Input;
use num_complex::Complex64;

// WE APPLY CONVENTION OF INDEX 0 HERE
pub fn apply_e_ij(m: &MetaVar, i: usize, j: usize) -> Option<(MetaVar, Complex64)> {
    if i >= m.len() || j >= m.len() { return None; }
    
    if m[j] == 0 { return None; }
    let mut m_copy = m.clone();
    let coeff = Complex64::new( m[i].into(),0.);
    if i == j { return Some((m_copy, coeff)); }

    m_copy[j] -= 1;
    m_copy[i] += 1;
    Some((m_copy, coeff+1.0))
}

fn apply_x_i_metavar(m: &MetaVar, vi: &Vec<usize>) -> Option<SLP> {
    unimplemented!()
}

use std::collections::HashMap;

pub fn apply_x_i(program: &SLP, vi: &Vec<usize>) -> Option<SLP> {
    if program.len() == 0 || vi.len() == 0 { return None;}

    let mut slp:SLP = Vec::new();
    let mut refs: HashMap<usize,usize> = HashMap::new();
    let l_line = program.last()?;


    for (i, line) in program.iter().enumerate() {
        match line {
            SLPLine::Input(m) => {
                slp.append( &mut apply_x_i_metavar(m, vi)? );
                refs.insert( i+1, slp.len());
            },
            SLPLine::Compound((s1, Plus, s2)) => {
                let mut s1 = *s1;
                let mut s2 = *s2;
                if let L(n) = s1 && refs.contains_key(&(n as usize)) {
                    s1 = L( refs[&(n as usize)] as u32);
                }
                if let L(n) = s2 && refs.contains_key(&(n as usize)) {
                    s2 = L( refs[&(n as usize)] as u32 );
                }
                slp.push( SLPLine::Compound((s1, Plus, s2)) );
                refs.insert( i+1, slp.len());
            }
            SLPLine::Compound((s1,Mult,s2)) => {}
        }
    }

    None

}
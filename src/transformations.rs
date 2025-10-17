use crate::straight_line_program::MetaVar;
use crate::straight_line_program::LineVar;
use crate::straight_line_program::Operation;
use crate::straight_line_program::SLPLine;
use crate::straight_line_program::SLPVar;
use crate::straight_line_program::SLP;

use num_complex::Complex64;

use std::collections::HashMap;
use std::collections::HashSet;

type TMatrix = Vec<u32>;
type LineMapper = HashMap<(TMatrix, LineVar), LineVar>;
type LinesToParse = HashMap<usize, HashSet<TMatrix>>;
type MetaVarSet = HashSet<MetaVar>;
type MetaVarLineMap = HashMap<MetaVar,usize>;


#[derive(Debug)]
enum PartialSLPVar {
    LineToTranslate( (TMatrix, LineVar) ),
    MetaVarRef( MetaVar ),
    Coeff( Complex64 ),
    LineInProgram( LineVar ),
    Zero
}
#[derive(Debug)]
enum PartialSLPLine {
    Zero,
    Input(MetaVar),
    Compound((PartialSLPVar, Operation, PartialSLPVar)),
}
type PartialSLP = Vec<PartialSLPLine>;


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

fn apply_x_i_input(m: &MetaVar, vi: &TMatrix, metavar_refs: &mut MetaVarSet) -> Option<(MetaVar, Complex64)> {
    // apply each of the bases for R^(k*k) one after another a number of times according to vi (ie. vis should have size k^2)
    // the bases need to be applied in lexigraphical ordering, an ordering not mentioned in the paper!
    // ORDERING ASSUMPTION: 
    //  will enumerate E_ij by i THEN j ie (0,0), (0,1), .., (0.k), (1,0), (1,1), ...
    // MULTIPLICATION ASSUMPTION:
    //  operation acts like differentiation by x_j multiplied by x_i => previous coefficients just combine:
    //  E_01 . 4*c<0,1,1> = 4 * E_01 . c<0,1,1> = 4 * c<1,0,1>

    let mut coeff = Complex64::ONE;
    let mut meta = m.clone();
    for i in 0..m.len() {
        for j in 0..m.len() {
            for _ in 0..vi[i*m.len() + j] {
                let (m,c) = apply_e_ij(&meta, i, j)?;
                meta = m;
                coeff *= c;
            }
        }
    }
    if coeff == Complex64::ZERO { return None;}
    metavar_refs.insert(meta.clone());
    Some((meta, coeff))
}

fn transform_x_i_input(m: &MetaVar, vi: &TMatrix, metavar_refs: &mut MetaVarSet, partial_program: &mut PartialSLP) {
    use PartialSLPLine::{Compound, Zero};
    use PartialSLPVar::{Coeff, MetaVarRef};
    use Operation::Mult;
    if let Some((meta, coeff)) = apply_x_i_input(m, vi, metavar_refs) {
        partial_program.push( Compound(( Coeff(coeff), Mult, MetaVarRef(meta) )) );
    } else {
        partial_program.push( Zero );
    }
}

fn add_matrix_to_lines_to_parse(line: LineVar, vi: TMatrix, lines_to_parse: &mut LinesToParse) {
    if let Some(hm) = lines_to_parse.get_mut(&line) {
        hm.insert( vi );
    } else {
        lines_to_parse.insert(line, HashSet::from([vi]) );
    }
}

fn transform_x_i_plus(s1: &SLPVar, s2: &SLPVar, vi: &TMatrix, partial_program: &mut PartialSLP, lines_to_parse: &mut LinesToParse ) {
    use PartialSLPLine::{Compound, Zero};
    use PartialSLPVar::{LineToTranslate, Coeff};
    use Operation::Plus;
    // We transform X^i . (a + b) = 0
    //              X^i . (Ln + b) = X^i . Ln
    //              X^i . (a + Ln) = X^i . Ln
    //              X^i . (Ln + Lm) = X^i . Ln + X^i . Lm

    let vi_is_zero =  vi.iter().all(|&n| n==0);

    if let SLPVar::F(c1) = *s1 && let SLPVar::F(c2) = *s2 {
        if vi_is_zero {
            partial_program.push( Compound(( Coeff(c1), Plus, Coeff(c2) )) );
        } else {
            partial_program.push( Zero );
        }
    }
    // TODO : if time allows, optimise this cloning of vectors away, could be very costly
    if let SLPVar::L(n) = *s1 && let SLPVar::F(c2) = *s2 {
        if vi_is_zero {
            partial_program.push( Compound(( Coeff(c2), Plus, LineToTranslate((vi.clone(), n)) )) );
        } else {
            partial_program.push( Compound(( PartialSLPVar::Zero, Plus, LineToTranslate((vi.clone(), n)) )) );
        }
        add_matrix_to_lines_to_parse(n, vi.clone(), lines_to_parse);
    }
    if let SLPVar::L(n) = *s2 && let SLPVar::F(c2) = *s1 {
        if vi_is_zero {
            partial_program.push( Compound(( Coeff(c2), Plus, LineToTranslate((vi.clone(), n)) )) );
        } else {
            partial_program.push( Compound(( PartialSLPVar::Zero, Plus, LineToTranslate((vi.clone(), n)) )) );
        }
        add_matrix_to_lines_to_parse(n, vi.clone(), lines_to_parse);
    }
    if let SLPVar::L(n1) = *s1 && let SLPVar::L(n2) = *s2 {
        partial_program.push( Compound(( LineToTranslate((vi.clone(), n1)), Plus, LineToTranslate((vi.clone(), n2)) )) );
        add_matrix_to_lines_to_parse(n1, vi.clone(), lines_to_parse);
        add_matrix_to_lines_to_parse(n2, vi.clone(), lines_to_parse);
    } 
}

fn transform_x_i_mult(s1: &SLPVar, s2: &SLPVar, vi: &TMatrix, partial_program: &mut PartialSLP, lines_to_parse: &mut LinesToParse ) {
    unimplemented!()
}

pub fn transform_x_i_program(init_program: &SLP, vi: &TMatrix) -> Option<SLP> {
    use SLPLine::*;
    use Operation::*;
    let mut metavars_used: MetaVarSet = MetaVarSet::new();
    let mut line_ref: LineMapper = LineMapper::new();
    let mut p_program: PartialSLP = PartialSLP::new();
    
    let xi_set = HashSet::from([ vi.clone()]);
    let mut line_parser: LinesToParse = HashMap::from([ (init_program.len()-1, xi_set ) ]);

    for (lno, line) in init_program.iter().enumerate().rev() {
        if let None = line_parser.get(&lno) {
            continue;
        }
        let transforms = line_parser.get(&lno).unwrap().clone();
        for tm in  transforms {
            line_ref.insert( (tm.clone(), lno), p_program.len() );
            match line {
                Input(m) => transform_x_i_input(m, &tm, &mut metavars_used, &mut p_program),
                Compound((s1, Plus, s2)) => transform_x_i_plus(s1, s2, &tm, &mut p_program, &mut line_parser),
                Compound((s1, Mult, s2)) => transform_x_i_mult(s1, s2, &tm, &mut p_program, &mut line_parser),
            }
        }
    }

    println!("metavars_used: {metavars_used:?}");
    println!("line_ref: {line_ref:?}");
    println!("line_parser: {line_parser:?}");
    println!("p_program: {p_program:?}");


    println!("\n STEP 2");
    line_ref.values_mut().map(|v| p_program.len() - v - 1 + metavars_used.len());
    let mut metavar_mapper = MetaVarLineMap::new();
    for (mi,m) in metavars_used.iter().enumerate() {
        metavar_mapper.insert(m.clone(), metavars_used.len()-1-mi);
        p_program.push(PartialSLPLine::Input(m.to_vec()));
    }
    p_program.reverse();

    use SLPVar::*;

    // TODO: Implement a Zero removing algorithm
    let mut program: SLP = SLP::new();
    for line in p_program {
        match line {
            // FIXME: temporary fix of setting zero values to zero + zero
            PartialSLPLine::Zero => program.push(Compound(( F(Complex64::ZERO), Plus, F(Complex64::ZERO) ))),
            PartialSLPLine::Input(m) => program.push( Input(m) ),
            PartialSLPLine::Compound((s1,op,s2)) => {

            }
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn apply_x_i_input_test() {
        // Hand calculated X^i . c<1,3,2>: X^i expressed as i = (0,0,1,0,2,0,0,0,0) in R^9 on c<1,3,2> should be 18*c<2,3,1>
        let expected = Some((vec![2,3,1], Complex64::ONE.scale(18.0) ));
        let mut hs = MetaVarSet::new();
        let result = apply_x_i_input(&vec![1,3,2], &vec![0,0,1,0,2,0,0,0,0],&mut hs);
        println!("metavars: {hs:?}");
        assert_eq!(expected, result);
    }

    #[test]
    fn project_x_i_plus_sanity() {
        use SLPLine::*;
        use SLPVar::*;
        use Operation::*;
        let program: SLP = vec![ Input(vec![1,3,2]), Input(vec![2,3,1]), Compound(( L(0), Plus, L(1) ))];

    }

    // TODO : add apply_x_i_input test for when vi is 0 vector
}
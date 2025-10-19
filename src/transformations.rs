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
type LineMapper = HashMap<LineVar, HashMap<TMatrix, LineVar>>;
type LinesToParse = HashMap<usize, HashSet<TMatrix>>;
type MetaVarSet = HashSet<MetaVar>;
type MetaVarLineMap = HashMap<MetaVar,usize>;


#[derive(Debug)]
enum PartialSLPVar {
    LineToTranslate( (TMatrix, LineVar) ),
    MetaVarRef( MetaVar ),
    Coeff( Complex64 ),
    LineInProgram( LineVar ), // for the multiplication case
}
#[derive(Debug)]
enum PartialSLPLine {
    Zero,
    Input(MetaVar),
    Compound((PartialSLPVar, Operation, PartialSLPVar)),
    Jump( (TMatrix, LineVar) ),
}
type PartialSLP = Vec<PartialSLPLine>;

fn stringify_partialslpvar(var: &PartialSLPVar) -> String {
    use PartialSLPVar::*;
    match var {
        LineToTranslate((tm,lno)) => format!("{tm:?}.L{lno}"),
        MetaVarRef(m) => format!("CR<{}>", m.into_iter().map(|i| i.to_string()).collect::<Vec<String>>().join(",")),
        Coeff(c) => c.to_string(),
        LineInProgram(lno) => format!("PL{lno}"),
    }
}
fn stringify_partialslpline(line: &PartialSLPLine) -> String {
    use PartialSLPLine::*;
    match line {
        Zero => "0".to_string(),
        Input(m) => format!("C<{}>",m.into_iter().map(|i| i.to_string()).collect::<Vec<String>>().join(",")),
        Compound((s1,Operation::Plus, s2)) => format!("{} + {}", stringify_partialslpvar(s1), stringify_partialslpvar(s2)),
        Compound((s1,Operation::Mult, s2)) => format!("{} * {}", stringify_partialslpvar(s1), stringify_partialslpvar(s2)),
        Jump((tm,lno)) => format!("{tm:?}.L{lno}"),
    }
}
pub fn stringify_partialslp(p_slp: &PartialSLP) -> String {
    p_slp.iter().enumerate().map(|(i,l)| "L".to_string() + i.to_string().as_str() + ": " + stringify_partialslpline(l).as_str() )
        .collect::<Vec<String>>().join("\n")
}


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
    use PartialSLPLine::{Compound, Zero, Jump};
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
            partial_program.push( Jump( (vi.clone(), n) ) ); 
        }
        add_matrix_to_lines_to_parse(n, vi.clone(), lines_to_parse);
    }
    if let SLPVar::L(n) = *s2 && let SLPVar::F(c2) = *s1 {
        if vi_is_zero {
            partial_program.push( Compound(( Coeff(c2), Plus, LineToTranslate((vi.clone(), n)) )) );
        } else {
            partial_program.push( Jump( (vi.clone(), n) ) ); 
        }
        add_matrix_to_lines_to_parse(n, vi.clone(), lines_to_parse);
    }
    if let SLPVar::L(n1) = *s1 && let SLPVar::L(n2) = *s2 {
        partial_program.push( Compound(( LineToTranslate((vi.clone(), n1)), Plus, LineToTranslate((vi.clone(), n2)) )) );
        add_matrix_to_lines_to_parse(n1, vi.clone(), lines_to_parse);
        add_matrix_to_lines_to_parse(n2, vi.clone(), lines_to_parse);
    } 
}

use std::cmp::min;
use std::fmt::format;
use std::iter::zip;

fn get_all_subvecs_i(v: &[u32], i: u32) -> Vec<TMatrix> {
    if i == 0 { return vec![Vec::from(v)]; }
    if i == v.iter().sum() { return vec![vec![0; v.len()]]; }
    if i > v.iter().sum() { return vec![]; }
    if v.len() == 0 { return vec![]; }
    if v.len() == 1 { return vec![vec![v[0]-i]]; }

    let mut subvecs_i = Vec::new();
    for j in 0..=min(i,v[0]) {
        let subvecs_i_j = get_all_subvecs_i(&v[1..], i-j);
        let mut subvecs_i_j: Vec<Vec<u32>> = subvecs_i_j.into_iter().map( |sv| [vec![v[0]-j], sv].concat() ).collect();
        subvecs_i.append(&mut subvecs_i_j);

    }
    subvecs_i
}
fn get_all_subvecs(v: &TMatrix) -> Vec<TMatrix> {
    (0..=v.iter().sum()).map(|i| get_all_subvecs_i(&v[..], i)).collect::<Vec<_>>().concat()
}

fn transform_x_i_mult(s1: &SLPVar, s2: &SLPVar, vi: &TMatrix, partial_program: &mut PartialSLP, lines_to_parse: &mut LinesToParse, lines_to_flip: &mut HashSet<usize> ) {
    use SLPVar::*;
    use PartialSLPVar::*;
    // Multiplication has complications:
    //  If both variables are coefficients then only result is Zero
    //  If one variable is a coefficient ( = Ln * c ) then return vi . Ln
    //  If both variables then nasty
    match (s1,s2) {
        (&F(_),&F(_)) => partial_program.push( PartialSLPLine::Zero ),
        (&F(_),&L(n)) => {
            partial_program.push( PartialSLPLine::Jump((vi.clone(), n)) );
            add_matrix_to_lines_to_parse(n, vi.clone(), lines_to_parse);
        },
        (&L(n),&F(_)) => {
            partial_program.push( PartialSLPLine::Jump((vi.clone(), n)) );
            add_matrix_to_lines_to_parse(n, vi.clone(), lines_to_parse);
        },
        (&L(n1),&L(n2)) => {
            let subvecs = get_all_subvecs(vi);
            println!("subvecs: {subvecs:?}");
            for (i,v_small) in get_all_subvecs(vi).into_iter().enumerate().rev() {

                let v_small_minus:TMatrix = zip(vi,&v_small).map(|(i,j)| i-j).collect();
                add_matrix_to_lines_to_parse(n1, v_small.clone(), lines_to_parse);
                add_matrix_to_lines_to_parse(n2, v_small_minus.clone(), lines_to_parse);
                if i != 0 {
                    partial_program.push( PartialSLPLine::Compound(( LineInProgram( partial_program.len()+1 ), Operation::Plus, LineInProgram( partial_program.len()+2 ))) );
                    lines_to_flip.insert(partial_program.len()-1);
                }
                partial_program.push( PartialSLPLine::Compound(( LineToTranslate((v_small,n1)), Operation::Mult, LineToTranslate((v_small_minus,n2)) )));
            }
        }

    }
}

fn insert_into_line_ref(line_ref: &mut LineMapper, lno: usize, tm: TMatrix, val: usize) {
    if let Some(tm_map) = line_ref.get_mut(&lno) {
        tm_map.insert( tm, val);    
    } else {
        let tm_map: HashMap<TMatrix, usize> = HashMap::from([(tm,val)]);
        line_ref.insert(lno, tm_map);
    }
}

pub fn transform_x_i_program(init_program: &SLP, vi: &TMatrix) -> Option<SLP> {
    use SLPLine::*;
    use Operation::*;
    let mut metavars_used: MetaVarSet = MetaVarSet::new();
    let mut line_ref: LineMapper = LineMapper::new();
    let mut p_program: PartialSLP = PartialSLP::new();
    
    let xi_set = HashSet::from([ vi.clone()]);
    let mut line_parser: LinesToParse = HashMap::from([ (init_program.len()-1, xi_set ) ]);

    let mut program_lines_to_flip: HashSet<usize> = HashSet::new();

    for (lno, line) in init_program.iter().enumerate().rev() {
        if let None = line_parser.get(&lno) {
            continue;
        }
        let transforms = line_parser.get(&lno).unwrap().clone();
        for tm in  transforms {
            insert_into_line_ref(&mut line_ref, lno, tm.clone(), p_program.len() );
            match line {
                Input(m) => transform_x_i_input(m, &tm, &mut metavars_used, &mut p_program),
                Compound((s1, Plus, s2)) => transform_x_i_plus(s1, s2, &tm, &mut p_program, &mut line_parser),
                Compound((s1, Mult, s2)) => transform_x_i_mult(s1, s2, &tm, &mut p_program, &mut line_parser, &mut program_lines_to_flip),
            }
        }
    }

    println!("metavars_used: {metavars_used:?}");
    println!("line_ref: {line_ref:?}");
    println!("line_parser: {line_parser:?}");
    println!("p_program: {p_program:?}");


    // Reverse all line reference values so they now point to where they will be in the list once reversed
    let flip_line_ref = |prog_len: usize, mvar_len: usize, index: LineVar| { prog_len - index - 1 + mvar_len };
    line_ref.values_mut().for_each(|m| {
        m.values_mut().for_each( |v| {*v = flip_line_ref(p_program.len(), metavars_used.len(), *v)} );
    });
    // Reverse all line references in the partial program also
    let p_program_len = p_program.len();
    program_lines_to_flip.iter().for_each(|lno| {
        use PartialSLPLine::*;
        use PartialSLPVar::*;
        if let Compound(tp) = &mut p_program[*lno] && let (LineInProgram(n1), _, LineInProgram(n2)) = tp {
            tp.0 = LineInProgram( flip_line_ref(p_program_len, metavars_used.len(), *n1));
            tp.2 = LineInProgram( flip_line_ref(p_program_len, metavars_used.len(), *n2));
        } else if let Compound(tp) = &mut p_program[*lno] && let (LineInProgram(n1), _, _) = tp {
            tp.0 = LineInProgram( flip_line_ref(p_program_len, metavars_used.len(), *n1));
        } else if let Compound(tp) = &mut p_program[*lno] && let (_, _, LineInProgram(n2)) = tp {
            tp.2 = LineInProgram( flip_line_ref(p_program_len, metavars_used.len(), *n2));
        }
    });

    // Chuck all the metavars onto the end of the program and add references to them in a map structure so references can be utilised
    let mut metavar_mapper = MetaVarLineMap::new();
    for (mi,m) in metavars_used.iter().enumerate() {
        metavar_mapper.insert(m.clone(), metavars_used.len()-1-mi);
        p_program.push(PartialSLPLine::Input(m.to_vec()));
    }
    p_program.reverse();

    // Section 2
    // 1. Translate all line references in the partial program to be in the right place
    // 2. Translate all (Matrix, LineNumber) references to be the correct line in the program
    // 3. Translate all references to metavariables to be the correct line in the program
    // 4. Apply a zero reduction algorithm to clear the amount of zeros that will probably clog up the output
    // 5. Convert the partial program (now only using inputs, zeros and compounds using coeffs and line references) into a program

    // step 1 (ignored until multiplication implemented)

    // step 2 & 3
    for line in p_program.iter_mut() {
        if let PartialSLPLine::Compound((s1, op, s2)) = line {
            use PartialSLPVar::*;
            let s1: PartialSLPVar = match s1 {
                LineInProgram(n) => LineInProgram(*n),
                Coeff(c) => Coeff(*c),
                LineToTranslate((tm, n)) => LineInProgram(line_ref[n][tm]),
                MetaVarRef(m) => LineInProgram( metavar_mapper[m] ),
            };
            let s2: PartialSLPVar = match s2 {
                LineInProgram(n) => LineInProgram(*n),
                Coeff(c) => Coeff(*c),
                LineToTranslate((tm, n)) => LineInProgram(line_ref[n][tm]),
                MetaVarRef(m) => LineInProgram( metavar_mapper[m] ),
            };
            *line = PartialSLPLine::Compound((s1, *op, s2));
        }
    }

    
    println!("Modified partial program:\n{}", stringify_partialslp(&p_program));


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
    fn get_all_subvecs_i_test() {
        let v = vec![1,2,3];
        let i = 2;
        let expected = HashSet::from( [ vec![0,1,3], vec![0,2,2], vec![1,0,3], vec![1,1,2], vec![1,2,1] ]) ;
        let result = get_all_subvecs_i(&v, i);
        let mut hash_res: HashSet<Vec<u32>> = HashSet::new();
        result.iter().for_each(|v| { hash_res.insert(v.clone()); });
        println!("got result: {result:?}");
        assert_eq!(expected, hash_res);
    }

    #[test]
    fn transform_x_i_mult_test1() {
        use SLPVar::L;
        let s1 = L(101);
        let s2 = L(202);
        let vi = vec![0,0,1,0,0,0,0,0,0];
        let mut p_program: PartialSLP = Vec::new();
        let mut lines_to_parse: LinesToParse = HashMap::new();
        let mut lines_to_flip: HashSet<usize> = HashSet::new();

        transform_x_i_mult(&s1, &s2, &vi, &mut p_program, &mut lines_to_parse, &mut lines_to_flip);
        println!("Partial Program: {p_program:?}");
        println!("lines to parse: {lines_to_parse:?}");
        println!("lines to flip: {lines_to_flip:?}");
        assert_eq!(1,1);
    }
}
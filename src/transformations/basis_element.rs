//! transformations provides the operations to apply the actions to a metapolynomial
//! Currently this only includes action by Eij, a kxk matrix where E(k,l) <- Indicator(k=i and l=j)
//! This element has a simpler action on a SLP and is all that is initially needed to program projections
//! If a constructive method to define a projection by a single element of GL_{k} is found, then the more general action may be implemented
#![allow(dead_code)]

use crate::straight_line_program::*;

use std::collections::HashMap;


/*
    To apply Eij action onto an SLP, the method works recursively backwards through the SLP
    Eij . m = m' ( or 0 ), declare 0 as 0+0 (as the SLP does not allow numbers to be a line by themselves) 
    Eij . (A+B) = (Eij . A) + (Eij . B)
    Eij . (A*B) = (Eij . A)*B + A*(Eij . B)

    Hence we construct the new SLP backwards and allow for references to operations not yet applied, like (Eij . L3)
    Going backwards through the SLP, a line may need to be included as is, or included transformed by Eij so store a dictionary <usize, <AsIs | Transformed | Both>>
     to see what needs to occur for each line
    Upon transforming a line, the result is stored in a PartialSLP, a structure designed to accomodate references to transformations not yet applied
     and the index of the transformed line is stored in another dictionary structure <(usize,bool), usize> where bool is whether this was the Transformed line or kept as it was
     to allow the later step to replace a transformation to be done with a line number reference to that new line.
    As the partial slp is built backwards and includes indices to itself, we choose to index as if it is the right way around and flip all LIP (Line references Inside Partial)
     at the end, as the final size won't be determined until the end

     Choices:
      - Are duplicates of metavars optimised away (before and/or after)
      - Are duplicates of lines optimised away (before and/or after)

*/

#[derive(Debug, PartialEq, Copy, Clone)]
enum Transforms {
    AsIs,
    Transformed,
    Both,
}

/// LinesToTransform holds dict of lines in original SLP and whether they are needed in transformed SLP, with or without transformation
type LinesToTransform = HashMap<usize, Transforms>;
/// LineMap holds dict of lines and transformation applied from original SLP and location in new transformed partial SLP
type LineMap = HashMap<(usize, bool), usize>;

/// PartialSLPVar is a structure mimicking SLPVar but in an unfinished format allowing it to be converted to a valid SLPVar later
pub enum PartialSLPVar<T> {
    LIP(usize), // LineInPartial is index of an earlier line in the transformed partial slp
    LTT((usize, bool)), // LineToTranslate, bool represents whether it should be transformed by eij or not
    C(T),
    Zero,
}
/// PartialSLPLine is a structure mimicking SLPLine but in an unfinished format allowing it to be converted to it later
pub enum PartialSLPLine<T> {
    Input(MetaVar),
    Compound((PartialSLPVar<T>, Operation, PartialSLPVar<T>)),
}
/// PartialSLPis a structure mimicking SLPbut in an unfinished format allowing it to be converted to it later
pub type PartialSLP<T> = Vec<PartialSLPLine<T>>;

use std::fmt::Display;
/// stringify_partialslpvar returns string presentation of partialslpvar for debugging
pub fn stringify_partialslpvar<T: Display>(var: &PartialSLPVar<T>) -> String {
    use PartialSLPVar::*;
    match var {
        LTT((lno, transforming)) => format!("LTT({lno},T?{transforming})"),
        C(c) => c.to_string(),
        LIP(lno) => format!("LIP{lno}"),
        Zero => format!("0"),
    }
}
/// stringify_partialslpline returns string presentation of partialslpline for debugging
pub fn stringify_partialslpline<T: Display>(line: &PartialSLPLine<T>) -> String {
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
pub fn stringify_partialslp<T: Display>(p_slp: &PartialSLP<T>) -> String {
    p_slp
        .iter()
        .enumerate()
        .map(|(i, l)| {
            "L".to_string() + i.to_string().as_str() + ": " + stringify_partialslpline(l).as_str()
        })
        .collect::<Vec<String>>()
        .join("\n")
}

/// apply_eij_on_metavar applies the group action of Eij onto a metavariable, as defined by Claim 4.2 of Algebraic metacomplexity and representation theory
pub fn apply_eij_on_metavar(m: &MetaVar, eij: (usize, usize)) -> Option<(u32, MetaVar)> {
    let (i, j) = eij;
    if i >= m.len() || j >= m.len() {
        // out of bounds for the metavariable => this is undefined behaviour / will result in 0
        None
    } else if m[j] == 0 {
        // always results in 0
        None
    } else if i != j {
        let mut m_new = m.clone();
        m_new[i] += 1;
        m_new[j] -= 1;
        Some((m[i] + 1, m_new))
    } else {
        Some((m[i], m.clone()))
    }
}

/// apply_eij_on_input_line computes the transformation of an input line (=C<1,2,3>) under eij
/// This involves updating the partial structures used to build the final slp: p_slp and line_map
/// m is the metavariable in the Input line, as a reference
/// lno is the line number of this Input line in the original SLP, needed to give the key for line_map
/// eij is the transformation being applied to the slp
/// p_slp is the partial slp structure being built to be converted into a slp, it is being built backwards
/// line_map is a dictionary mapping original line numbers and to their current index in the partial slp structure
/// u32_to_c is a conversion function to allow the coefficients from the metavar to be included in the program, a sensible definition should have u32_to_c(0) == T::default()
fn apply_eij_on_input_line<T, F: Fn(u32) -> T>(
    m: &MetaVar,
    lno: usize,
    eij: (usize, usize),
    p_slp: &mut PartialSLP<T>,
    line_map: &mut LineMap,
    u32_to_c: F,
) {
    use Operation::*;
    use PartialSLPLine::*;
    use PartialSLPVar::*;
    if let Some((coeff, m_new)) = apply_eij_on_metavar(m, eij) {
        // result is of form coeff * m_new, expressing in SLP is =m_new \n =Coeff*L(before)
        // partial_slp is constructed backwards meaning we add as follows: =Coeff*L(current slp len+1), =m_new
        //  the +1 is necessary as upon pushing, current slp len refers to index of line to be added, we want to refer to the line afterwards
        let p_slp_len = p_slp.len();
        p_slp.push(Compound((C(u32_to_c(coeff)), Mult, LIP(p_slp_len + 1))));
        p_slp.push(Input(m_new));
        // add reference to the new entry point for Eij . m into the line mapper
        line_map.insert((lno, true), p_slp_len);
    } else {
        // result is 0
        // know that by rules of SLP, a line cannot refer to itself => in p_slp will never have variable of LIP(0) hence this is reserved for 0
        line_map.insert((lno, true), 0);
    }
}
/// add_to_line_to_trans adds a line to be computed later by the transformation of the slp
///  as it can be either included as is or transformed or both, these cases need to be addressed when changing the structure
///  this function modifies lines_to_trans to update the hashmap with the correct transformations needed for lno
/// lno is the line to be included (transformed or not)
/// transform is the current transformation of lno desired
/// lines_to_transform is the hashmap of the lines and their transformations needed
fn add_to_line_to_trans(lno: usize, transform: Transforms, lines_to_trans: &mut LinesToTransform) {
    if let Some(t) = lines_to_trans.get(&lno) {
        if t != &transform {
            lines_to_trans.insert(lno, Transforms::Both);
        } else {
            lines_to_trans.insert(lno, transform);
        }
    } else {
        lines_to_trans.insert(lno, transform);
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
fn apply_eij_on_addition<T>(
    v1: &SLPVar<T>,
    v2: &SLPVar<T>,
    lno: usize,
    p_slp: &mut PartialSLP<T>,
    lines_to_trans: &mut LinesToTransform,
    line_map: &mut LineMap,
) {
    use Operation::*;
    use PartialSLPLine::*;
    use PartialSLPVar::{LTT, Zero};
    use SLPVar::*;
    match (v1, v2) {
        (L(n1), L(n2)) => {
            add_to_line_to_trans(*n1, Transforms::Transformed, lines_to_trans);
            add_to_line_to_trans(*n2, Transforms::Transformed, lines_to_trans);
            p_slp.push(Compound((LTT((*n1, true)), Plus, LTT((*n2, true)))));
        }
        (L(n), C(_)) | (C(_), L(n)) => {
            add_to_line_to_trans(*n, Transforms::Transformed, lines_to_trans);
            p_slp.push(Compound((LTT((*n, true)), Plus, Zero)))
        }
        (C(_), C(_)) => p_slp.push(Compound((Zero, Plus, Zero))),
    }
    line_map.insert((lno, true), p_slp.len() - 1);
}
/// apply_eij_on_product computes the transformation of a product line (=L2*0.5) under eij
/// This involves updating the partial structures used to build the final slp: p_slp and line_map
/// v1 the the first term in the line
/// v2 the the second term in the line
/// lno is the line number of this Input line in the original SLP, needed to give the key for line_map
/// p_slp is the partial slp structure being built to be converted into a slp, it is being built backwards
/// lines_to_trans is a dictionary of other lines that also need to be included into the partial slp
/// line_map is a dictionary mapping original line numbers and to their current index in the partial slp structure
fn apply_eij_on_product<T: Clone>(
    v1: &SLPVar<T>,
    v2: &SLPVar<T>,
    lno: usize,
    p_slp: &mut PartialSLP<T>,
    lines_to_trans: &mut LinesToTransform,
    line_map: &mut LineMap,
) {
    use Operation::*;
    use PartialSLPLine::*;
    use PartialSLPVar::C as PC;
    use PartialSLPVar::{LIP, LTT, Zero};
    use SLPVar::*;
    let p_len = p_slp.len();
    match (v1, v2) {
        (L(n1), L(n2)) => {
            add_to_line_to_trans(*n1, Transforms::Both, lines_to_trans);
            add_to_line_to_trans(*n2, Transforms::Both, lines_to_trans);
            p_slp.push(Compound((LIP(p_len + 1), Plus, LIP(p_len + 2))));
            p_slp.push(Compound((LTT((*n1, true)), Mult, LTT((*n2, false)))));
            p_slp.push(Compound((LTT((*n1, false)), Mult, LTT((*n2, true)))));
        }
        (L(n), C(c)) | (C(c), L(n)) => {
            add_to_line_to_trans(*n, Transforms::Transformed, lines_to_trans);
            p_slp.push(Compound((PC(c.clone()), Mult, LTT((*n, true)))));
        }
        (C(_), C(_)) => p_slp.push(Compound((Zero, Plus, Zero))),
    }
    line_map.insert((lno, true), p_len);
}
/// apply_eij_on_line computes the transformation of a slp line
/// This involves matching the line to the correct form to invoke the correct transformation function
/// line the the line being computed
/// lno is the line number of this Input line in the original SLP, needed to give the key for line_map
/// p_slp is the partial slp structure being built to be converted into a slp, it is being built backwards
/// lines_to_trans is a dictionary of other lines that also need to be included into the partial slp
/// line_map is a dictionary mapping original line numbers and to their current index in the partial slp structure
fn apply_eij_to_line<T, F>(
    line: &SLPLine<T>,
    lno: usize,
    eij: (usize, usize),
    p_slp: &mut PartialSLP<T>,
    lines_to_trans: &mut LinesToTransform,
    line_map: &mut LineMap,
    u32_to_c: F,
) where
    T: Clone,
    F: Fn(u32) -> T,
{
    use Operation::*;
    use SLPLine::*;
    match line {
        Input(m) => apply_eij_on_input_line(m, lno, eij, p_slp, line_map, u32_to_c),
        Compound((v1, Plus, v2)) => {
            apply_eij_on_addition(v1, v2, lno, p_slp, lines_to_trans, line_map)
        }
        Compound((v1, Mult, v2)) => {
            apply_eij_on_product(v1, v2, lno, p_slp, lines_to_trans, line_map)
        }
    }
}
/// include_line_from_slp includes a line from the SLP as is into the partial slp
/// This means copying the line across and ensuring other lines referenced will be copied across also
/// line the the line being computed
/// lno is the line number of this Input line in the original SLP, needed to give the key for line_map
/// p_slp is the partial slp structure being built to be converted into a slp, it is being built backwards
/// lines_to_trans is a dictionary of other lines that also need to be included into the partial slp
/// line_map is a dictionary mapping original line numbers and to their current index in the partial slp structure
fn include_line_from_slp<T: Clone>(
    line: &SLPLine<T>,
    lno: usize,
    p_slp: &mut PartialSLP<T>,
    lines_to_trans: &mut LinesToTransform,
    line_map: &mut LineMap,
) {
    use PartialSLPLine::Compound as PCompound;
    use PartialSLPLine::Input as PInput;
    use PartialSLPVar::C as PC;
    use PartialSLPVar::LTT;
    use SLPLine::*;
    use SLPVar::*;
    let p_len = p_slp.len();
    match line {
        Compound((L(n1), op, L(n2))) => {
            add_to_line_to_trans(*n1, Transforms::AsIs, lines_to_trans);
            add_to_line_to_trans(*n2, Transforms::AsIs, lines_to_trans);
            p_slp.push(PCompound((LTT((*n1, false)), *op, LTT((*n2, false)))));
        }
        Compound((L(n), op, C(c))) | Compound((C(c), op, L(n))) => {
            add_to_line_to_trans(*n, Transforms::AsIs, lines_to_trans);
            p_slp.push(PCompound((PC(c.clone()), *op, LTT((*n, false)))));
        }
        Compound((C(c1), op, C(c2))) => {
            p_slp.push(PCompound((PC(c1.clone()), *op, PC(c2.clone()))))
        }
        Input(m) => p_slp.push(PInput(m.clone())),
    }
    line_map.insert((lno, false), p_len);
}

/// apply_eij_on_program computes the transformation applied to the program
/// slp is the straight line program
/// eij is the transformation
/// u32_to_c is a function converting a u32 to the custom coefficient type T, ensure that u32_to_c(0) is sensible (ie. is 0 in your vector field)
/// returns result of slp or error message about failure
pub fn apply_eij_on_program<T, F>(
    slp: &SLP<T>,
    eij: (usize, usize),
    u32_to_c: F,
) -> Result<SLP<T>, String>
where
    T: Clone + Display,
    F: Fn(u32) -> T + Copy,
{
    use Operation::*;
    use PartialSLPLine::*;
    use PartialSLPVar::*;
    // Steps to convert SLP:
    //  0. Validity of inputs
    //  1. Initialisation of structures
    //  2. Calculation of each line
    //  3. Invert PartialSLP, LTT -> LIP,  correction of indices, conversion to SLP

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
    if eij.0 >= get_metavar_len(slp).unwrap() || eij.1 >= get_metavar_len(slp).unwrap() {
        return Err("Eij refers to a matrix larger than the size of the metavariables".into());
    }

    // 1. Initialisation
    let mut p_slp: PartialSLP<T> = vec![];
    let mut lines_to_trans: LinesToTransform = HashMap::new();
    let mut line_map: LineMap = HashMap::new();

    // include the final line as the first line to transform, to propagate backwards through the slp
    lines_to_trans.insert(slp.len() - 1, Transforms::Transformed);

    // 2. Calculation of each line
    for (lno, line) in slp.iter().enumerate().rev() {
        if let Some(&t) = lines_to_trans.get(&lno) {
            if t == Transforms::AsIs || t == Transforms::Both {
                include_line_from_slp(line, lno, &mut p_slp, &mut lines_to_trans, &mut line_map);
            }
            if t == Transforms::Transformed || t == Transforms::Both {
                apply_eij_to_line(
                    line,
                    lno,
                    eij,
                    &mut p_slp,
                    &mut lines_to_trans,
                    &mut line_map,
                    u32_to_c,
                );
            }
        }
    }

    // 3. Inversion and conversion
    p_slp.push(Compound((C(u32_to_c(0)), Plus, C(u32_to_c(0))))); // add zero to the beginning (once reversed)
    p_slp.reverse();

    // convert all LTT to LIP as at this point every line should have a transformation in the p_slp
    for line in p_slp.iter_mut() {
        match line {
            Compound((v1, _, v2)) => {
                if let LTT(p1) = v1 {
                    *v1 = LIP(*line_map.get(p1).ok_or("Couldn't find transformation")?);
                }
                if let LTT(p2) = v2 {
                    *v2 = LIP(*line_map.get(p2).ok_or("Couldn't find transformation")?);
                }
            }
            Input(_) => (),
        }
    }

    // convert all Zero to a reference to 0
    for line in p_slp.iter_mut() {
        if let Compound((v1, _, v2)) = line {
            if matches!(v1, Zero) {
                *v1 = LIP(0);
            }
            if matches!(v2, Zero) {
                *v2 = LIP(0);
            }
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

    // convert PartialSLP to SLP, should be case that all variables are LIP or C
    use SLPLine::Compound as SLPCompound;
    use SLPLine::Input as SLPInput;
    use SLPVar::C as SLPC;
    use SLPVar::L;
    let mut slp: SLP<T> = Vec::new();
    let pv_to_sv = |v: PartialSLPVar<T>| match v {
        LIP(n) => Ok(L(n)),
        C(c) => Ok(SLPC(c.clone())),
        LTT(_) => Err("Unconverted LTT found"),
        Zero => Err("Unconverted Zero found"),
    };
    for line in p_slp {
        slp.push(match line {
            Compound((v1, op, v2)) => SLPCompound((pv_to_sv(v1)?, op, pv_to_sv(v2)?)),
            Input(m) => SLPInput(m.clone()),
        });
    }
    Ok(slp)
}

/// apply_eijs_on_program computes the transformations applied to the program in REVERSE order (as function application is carried out that way)
/// slp is the straight line program
/// eijs is the vector of transformations, to be applied in order of last to first
/// u32_to_c is a function converting a u32 to the custom coefficient type T, ensure that u32_to_c(0) is sensible (ie. is 0 in your vector field)
/// returns result of slp or error message about first failure
pub fn apply_eijs_on_program<T, F>(
    slp: &SLP<T>,
    eijs: &Vec<(usize, usize)>,
    u32_to_c: F,
) -> Result<SLP<T>, String>
where
    T: Clone + Display,
    F: Fn(u32) -> T + Copy,
{
    eijs.iter()
        .enumerate()
        .rev()
        .try_fold(slp.clone(), |s, (i, &eij)| {
            apply_eij_on_program(&s, eij, u32_to_c).map_err(|err| {
                format!(
                    "Failure index {i}, Error: {err}, SLP: {}",
                    stringify_slp(&s)
                )
            })
        })
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::evaluation;
    use crate::parsing;
    use num_rational::Rational64;

    #[test]
    fn apply_x_i_test() {
        // let slp_str = "=C<1,2,3,4>\n=L0*L0";
        let slp_str = "=C<0,0,2>
=C<0,2,0>
=C<2,0,0>
=L0*L1
=L3*L2";
        let eij = (0, 0);
        let slp = parsing::slp_parser_rational(slp_str).expect("Should parse SLP");
        let slp_res = apply_eij_on_program(&slp, eij, |n| Rational64::from_integer(n.into()))
            .expect("Should apply Eij on SLP");

        println!(
            "Transformation: {eij:?} Parsed SLP:{}\n",
            stringify_slp(&slp)
        );
        println!("Result SLP:{}\n", stringify_slp(&slp_res));

        let dist = evaluation::RationalDistr {};
        let mut rng = rand::rng();
        let metavar_point = evaluation::get_random_val_for_metavars(&slp_res, &dist, &mut rng);

        let expected_slp_str = "=C<0,0,2>
=C<0,2,0>
=C<2,0,0>
=L0*L1
=L3*L2
=2*L4";
        let expected_slp =
            parsing::slp_parser_rational(expected_slp_str).expect("Should parse expected SLP");

        println!("Expected SLP:{}\n", stringify_slp(&expected_slp));
        println!("Evaluated at: {metavar_point:?}");

        assert!(evaluation::are_slps_similar(
            &metavar_point,
            &slp_res,
            &expected_slp
        ));
    }
}

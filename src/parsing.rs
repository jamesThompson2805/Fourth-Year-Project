//! parsing provides utilities for parsing SLP from strings
//! Format of a SLP string should be lines of either lists in form C<a,b,c,d> where a,b,c,d are non-neg integers representing metavariables
//!  or addition/multiplication expressions composed of coefficients or line references
//!  coefficients are of form determined by choice of slp_parser_..., see coeff_parser_... for formatting
//!  line references are of form Ln where n is non-neg integer, parser enforces that n must be less than the current line number (line numbers are 0-indexed)
//! There is no whitespace permitted, all lines start with = to aid readability in text file
//! eg.  =C<10,0,3> defines a metavariable
//! eg2. =L3+0.4 defines an addition operation, adding 0.4 to the result of line 3
//! eg3. =L4*L5 defines multiplication of the result of line 4 to line 5
//! The SLP is evaluated as the final line being the end result of the SLP, so care should be taken to ensure that you make no mistakes defining it
//! =C<0,1,2>
//! =L0*L0
//! =L1*L1
//! =5+5
//! The SLP above would be treated as the equation 5+5, and would ignore all previous lines
#![allow(dead_code)]
use super::straight_line_program::*;

use nom::bytes::complete::tag;
use nom::character::complete::digit1;
use nom::combinator::{map, map_res, value};
use nom::error::Error;
use nom::multi::separated_list1;
use nom::Parser;
use nom::IResult;
use nom::branch::alt;
use nom::number::complete::float;

use num_complex::Complex64;
use num_rational::Rational64;

/// Parses u32 decimal numbers
fn dec32_parser(input: &str) -> IResult<&str,u32> {
    map_res( digit1, |s: &str| s.parse::<u32>()).parse(input)
}
/// Parses usize decimal numbers
fn dec_parser(input: &str) -> IResult<&str,usize> {
    map_res( digit1, |s: &str| s.parse::<usize>()).parse(input)
}

/// Parses lists of u32 decimal numbers
fn metavar_parser(input: &str) -> IResult<&str, MetaVar> {
    map( (tag("C"),tag("<"), separated_list1(tag(","), dec32_parser), tag(">")),
     |c| c.2).parse(input)
}

/// Parses line references (Ln where n is a usize number) 
fn linevar_parser(input: &str) -> IResult<&str, usize> {
    map( (tag("L"), dec_parser),
     |c| c.1).parse(input)
}

/// Parses complex coefficients in format (a+bi) where a and b are f32, returns a complex number of f64s
fn coeff_parser_complex(input: &str) -> IResult<&str, Complex64> {
    map( (tag("("),float, tag("+"), float, tag("i)")), 
    |t| Complex64::new(t.1.into(),t.3.into())).parse(input)
}
/// Parses real coefficients into f64
fn coeff_parser_real(input: &str) -> IResult<&str, f64> {
    float.parse(input).map(|(s,f)| (s,f as f64))
}

/// Parses rational coefficients into a rational number format, parses from float representation or fraction form
fn coeff_parser_fraction(input: &str) -> IResult<&str, Rational64> {
    map( (tag("("),dec32_parser, tag("/"), dec32_parser,tag(")")), |t| Rational64::new(t.1.into(), t.3.into())).parse(input)
}

/// Parses rational coefficients into a rational number format, parses from float representation so may be inaccuracy or failure
fn coeff_parser_rational_from_float(input: &str) -> IResult<&str, Rational64> {
    map_res( coeff_parser_real, |f|
        Rational64::approximate_float(f).ok_or("Couldn't parse float as rational")
    ).parse(input)
}

/// Parses rational coefficients into a rational number format, parses from float representation or fraction form
fn coeff_parser_rational(input: &str) -> IResult<&str, Rational64> {
    alt((coeff_parser_fraction, coeff_parser_rational_from_float)).parse(input)
}

fn var_parser_complex(input: &str) -> IResult<&str, SLPVar<Complex64>> {
    alt((
        map( coeff_parser_complex, SLPVar::C), 
        map( linevar_parser, SLPVar::L)
    )).parse(input)
}
fn var_parser_real(input: &str) -> IResult<&str, SLPVar<f64>> {
    alt((
        map( coeff_parser_real, SLPVar::C), 
        map( linevar_parser, SLPVar::L)
    )).parse(input)
}
fn var_parser_rational(input: &str) -> IResult<&str, SLPVar<Rational64>> {
    alt((
        map( coeff_parser_rational, SLPVar::C), 
        map( linevar_parser, SLPVar::L)
    )).parse(input)
}

fn op_parser(input: &str) -> IResult<&str, Operation> {
    alt((
        value(Operation::Plus, tag("+")),
        value(Operation::Mult, tag("*"))
    )).parse(input)
}

fn line_parser_complex(input:&str) -> IResult<&str, SLPLine<Complex64>> {
    let (input,_) = tag("=")(input)?;
    alt((
        map((var_parser_complex, op_parser, var_parser_complex), |t|SLPLine::Compound(t)),
        map( metavar_parser, |v| SLPLine::Input(v)),
    )).parse(input)
}
fn line_parser_real(input:&str) -> IResult<&str, SLPLine<f64>> {
    let (input,_) = tag("=")(input)?;
    alt((
        map((var_parser_real, op_parser, var_parser_real), |t|SLPLine::Compound(t)),
        map( metavar_parser, |v| SLPLine::Input(v)),
    )).parse(input)
}
fn line_parser_rational(input:&str) -> IResult<&str, SLPLine<Rational64>> {
    let (input,_) = tag("=")(input)?;
    alt((
        map((var_parser_rational, op_parser, var_parser_rational), |t|SLPLine::Compound(t)),
        map( metavar_parser, |v| SLPLine::Input(v)),
    )).parse(input)
}

use nom::Err as Err2;

pub fn slp_parser_complex(input: &str) -> Result<SLP<Complex64>,String> {
    let mut program: SLP<Complex64> = Vec::new();
    for (l_no,line) in input.split('\n').enumerate() {
        let res = line_parser_complex(line);
        match res {
            Ok((_s,slp_line)) => {
                //if s != "" { return Err( format!("Error: line {l_no} : {line} not fully parsed, didn't parse {s}")); }
                if let SLPLine::Compound(slp_tuple) = &slp_line {
                    if let SLPVar::L(n) = slp_tuple.0 && n >= l_no.try_into().unwrap()  { return Err( format!("Error: line {l_no} : {line} refering to undeclared line : L{n}")); }
                    if let SLPVar::L(n) = slp_tuple.2 && n >= l_no.try_into().unwrap()  { return Err( format!("Error: line {l_no} : {line} refering to undeclared line : L{n}")); }
                }

                program.push(slp_line);
            },
            Err(Err2::Incomplete(_)) => return Err( format!("Error: line {l_no} : {line} is incomplete")),
            Err(Err2::Error(Error {input: _, code})) => return Err( format!("Error: line {l_no} : {line} trying to parse {code:?} ") ),
            Err(Err2::Failure(Error {input: _, code})) => return Err( format!("Error: line {l_no} : {line} trying to parse {code:?} ") ),
        }
    }
    Ok(program)
}
pub fn slp_parser_real(input: &str) -> Result<SLP<f64>,String> {
    let mut program: SLP<f64> = Vec::new();
    for (l_no,line) in input.split('\n').enumerate() {
        let res = line_parser_real(line);
        match res {
            Ok((_s,slp_line)) => {
                //if s != "" { return Err( format!("Error: line {l_no} : {line} not fully parsed, didn't parse {s}")); }
                if let SLPLine::Compound(slp_tuple) = &slp_line {
                    if let SLPVar::L(n) = slp_tuple.0 && n >= l_no.try_into().unwrap()  { return Err( format!("Error: line {l_no} : {line} refering to undeclared line : L{n}")); }
                    if let SLPVar::L(n) = slp_tuple.2 && n >= l_no.try_into().unwrap()  { return Err( format!("Error: line {l_no} : {line} refering to undeclared line : L{n}")); }
                }

                program.push(slp_line);
            },
            Err(Err2::Incomplete(_)) => return Err( format!("Error: line {l_no} : {line} is incomplete")),
            Err(Err2::Error(Error {input: _, code})) => return Err( format!("Error: line {l_no} : {line} trying to parse {code:?} ") ),
            Err(Err2::Failure(Error {input: _, code})) => return Err( format!("Error: line {l_no} : {line} trying to parse {code:?} ") ),
        }
    }
    Ok(program)
}
pub fn slp_parser_rational(input: &str) -> Result<SLP<Rational64>,String> {
    let mut program: SLP<Rational64> = Vec::new();
    for (l_no,line) in input.split('\n').enumerate() {
        let res = line_parser_rational(line);
        match res {
            Ok((_s,slp_line)) => {
                //if s != "" { return Err( format!("Error: line {l_no} : {line} not fully parsed, didn't parse {s}")); }
                if let SLPLine::Compound(slp_tuple) = &slp_line {
                    if let SLPVar::L(n) = slp_tuple.0 && n >= l_no.try_into().unwrap()  { return Err( format!("Error: line {l_no} : {line} refering to undeclared line : L{n}")); }
                    if let SLPVar::L(n) = slp_tuple.2 && n >= l_no.try_into().unwrap()  { return Err( format!("Error: line {l_no} : {line} refering to undeclared line : L{n}")); }
                }

                program.push(slp_line);
            },
            Err(Err2::Incomplete(_)) => return Err( format!("Error: line {l_no} : {line} is incomplete")),
            Err(Err2::Error(Error {input: _, code})) => return Err( format!("Error: line {l_no} : {line} trying to parse {code:?} ") ),
            Err(Err2::Failure(Error {input: _, code})) => return Err( format!("Error: line {l_no} : {line} trying to parse {code:?} ") ),
        }
    }
    Ok(program)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evaluation::are_rational_slp_similar;

    #[test]
    fn test_parsing_rationals() {
        let slp_str = "=C<0,0,1>
=L0*1.5";
        let slp_str2 = "=C<0,0,1>
=L0*(3/2)";
        let slp1 = slp_parser_rational( slp_str ).unwrap();
        let slp2 = slp_parser_rational( slp_str2 ).unwrap();
        assert!(are_rational_slp_similar(&slp1, &slp2).1);
    }
}
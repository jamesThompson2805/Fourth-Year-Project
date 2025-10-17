use crate::straight_line_program::MetaVar;
use crate::straight_line_program::LineVar;
use crate::straight_line_program::Operation;
use crate::straight_line_program::SLPLine;
use crate::straight_line_program::SLPVar;
use crate::straight_line_program::SLP;

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

fn dec32_parser(input: &str) -> IResult<&str,u32> {
    map_res( digit1, |s: &str| s.parse::<u32>()).parse(input)
}
fn dec_parser(input: &str) -> IResult<&str,usize> {
    map_res( digit1, |s: &str| s.parse::<usize>()).parse(input)
}

fn metavar_parser(input: &str) -> IResult<&str, MetaVar> {
    map( (tag("C"),tag("<"), separated_list1(tag(","), dec32_parser), tag(">")),
     |c| c.2).parse(input)
}

fn linevar_parser(input: &str) -> IResult<&str, LineVar> {
    map( (tag("L"), dec_parser),
     |c| c.1).parse(input)
}

fn coeff_parser(input: &str) -> IResult<&str, Complex64> {
    map( (tag("("),float, tag("+"), float, tag("i)")), 
    |t| Complex64::new(t.1.into(),t.3.into())).parse(input)
}


fn var_parser(input: &str) -> IResult<&str, SLPVar> {
    alt((
        map( coeff_parser, SLPVar::F), 
        map( linevar_parser, SLPVar::L)
    )).parse(input)
}

fn op_parser(input: &str) -> IResult<&str, Operation> {
    alt((
        value(Operation::Plus, tag("+")),
        value(Operation::Mult, tag("*"))
    )).parse(input)
}

fn line_parser(input:&str) -> IResult<&str, SLPLine> {
    let (input,_) = tag("=")(input)?;
    alt((
        map((var_parser, op_parser, var_parser), |t|SLPLine::Compound(t)),
        map( metavar_parser, |v| SLPLine::Input(v)),
    )).parse(input)
}

use nom::Err as Err2;

pub fn file_parser(input: &str) -> Result<SLP,String> {
    let mut program: SLP = Vec::new();
    for (l_no,line) in input.split('\n').enumerate() {
        let res = line_parser(line);
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
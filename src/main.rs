use crate::parsing::slp_parser_rational;
use crate::projections::apply_lambda_projection_to_slp;
use crate::straight_line_program::{generate_random_homogeneous_slp, reduce_slp};

mod evaluation;
mod parsing;
mod pbw_reduce;
mod projections;
mod straight_line_program;
mod transformations;

use straight_line_program::*;
use evaluation::*;
use parsing::*;

use num_rational::Rational64;
use num_bigint::BigInt;
use rand::{Rng, rng};

use std::fmt::{Display, Debug};
use std::fs;
use std::fs::File;
use std::io::Write;


fn save_slp_to_file<T: Display + Debug>(slp: &SLP<T>, out_name: String) {
    let path = out_name+".txt";
    let mut output = File::create(path).unwrap();
    write!(output, "{}", stringify_slp(&slp)).unwrap();
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    
    let mut coeff_rng = rng();

    let mut gen_coeff = || {
        let numer: i64 = coeff_rng.random_range(-10..10);
        let denom: i64 = coeff_rng.random_range(1..10);  // strictly positive to avoid division by zero
        Rational64::new(numer, denom)
    };
    let i64_to_c = |i| Rational64::new(i, 1);
    let bigi_to_c = |i:BigInt| Rational64::new(i.try_into().unwrap(), 1);

    let slp_zero_str = "=0+0";
    let slp_zero = slp_parser_rational(slp_zero_str)?;

    let mut found = false;
    while !found {
        let mut slp = generate_random_homogeneous_slp::<Rational64,_,_>(&mut gen_coeff, 11, 3, 4, &mut rng());
        reduce_slp(&mut slp, Rational64::ZERO, Rational64::ONE);

        let slp_res = apply_lambda_projection_to_slp::<_,_,16>(&slp, &vec![6,6,6,6], bigi_to_c)?;

        if !evaluation::are_rational_slp_similar(&slp_res, &slp_zero).1 {
            save_slp_to_file(&slp_res, "11_3_4_proj_found.txt".into());
            found = true;
        }
    }

    Ok(())
}

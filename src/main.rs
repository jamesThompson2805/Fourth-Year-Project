mod straight_line_program;
mod parsing;
mod transformations;

use std::fs;

use crate::{straight_line_program::stringify_slp, transformations::transform_x_i_program};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let contents = fs::read_to_string("input.txt")
        .expect("Should have been able to read the file");

    let res = parsing::file_parser(&contents);

   match res {
    Ok(slp) => {
        println!("gathered output: {slp:?}");
        println!();
        let maybe_slp = transform_x_i_program(&slp, &vec![3,0,0,3,0,0,3,0,0]);
        if let Some(slp) = maybe_slp {
            println!("Success, received SLP:\n{}", stringify_slp(&slp));
        }
    },

    Err(s) => println!("ERROR IN PARSING: \n{s}"),
   } 

    Ok(())
}

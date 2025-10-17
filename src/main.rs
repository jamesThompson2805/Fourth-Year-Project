mod straight_line_program;
mod parsing;
mod transformations;

use std::fs;

use crate::transformations::transform_x_i_program;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let contents = fs::read_to_string("input.txt")
        .expect("Should have been able to read the file");

    let res = parsing::file_parser(&contents);

   match res {
    Ok(slp) => {
        println!("gathered output: {slp:?}");
        println!();
        transform_x_i_program(&slp, &vec![0,0,1,0,0,0,0,0,0]);
    },

    Err(s) => println!("ERROR IN PARSING: \n{s}"),
   } 

    Ok(())
}

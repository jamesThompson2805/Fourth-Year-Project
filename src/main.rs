mod straight_line_program;
mod parsing;
mod transformations;
mod projections;

use std::fs;
use std::fs::File;
use std::io::{Error, Write};

use crate::{straight_line_program::stringify_slp, transformations::transform_x_i_program};

fn save_input_to_output(vi: &Vec::<u32>,out_name: String) -> Result<(), Box<dyn std::error::Error>> {
    let contents = fs::read_to_string("input.txt")
        .expect("Should have been able to read the file");

    let res = parsing::file_parser(&contents);

   match res {
    Ok(slp) => {
        println!("gathered output: {slp:?}");
        println!();
        let maybe_slp = transform_x_i_program(&slp, vi);
        if let Some(res_slp) = maybe_slp {
            println!("Success, received SLP");
            let path = out_name+".txt";
            let mut output = File::create(path)?;
            write!(output, "===== INPUT SLP =====\n{}\n===== TRANSFORM =====\n{:?}\n===== OUTPUT SLP =====\n{}", stringify_slp(&slp), vi, stringify_slp(&res_slp))?;
        }
    },

    Err(s) => println!("ERROR IN PARSING: \n{s}"),
   } 
   Ok(())
}

fn print_input_to_output(vi: &Vec::<u32>) -> Result<(), Box<dyn std::error::Error>> {
    let contents = fs::read_to_string("input.txt")
        .expect("Should have been able to read the file");

    let res = parsing::file_parser(&contents);

   match res {
    Ok(slp) => {
        println!("gathered output: {slp:?}");
        println!();
        let maybe_slp = transform_x_i_program(&slp, vi);
        if let Some(slp) = maybe_slp {
            println!("Success, received SLP:\n{}", stringify_slp(&slp));
        }
    },

    Err(s) => println!("ERROR IN PARSING: \n{s}"),
   } 

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let vi = vec![3,0,0,3,0,0,3,0,0];
    save_input_to_output(&vi,"example_5_4".to_string())
}

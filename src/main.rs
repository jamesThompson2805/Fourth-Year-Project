mod straight_line_program;
mod parsing;
mod projection;

use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let contents = fs::read_to_string("input.txt")
        .expect("Should have been able to read the file");

    let res = parsing::file_parser(&contents);

   match res {
    Ok(slp) => println!("gathered output: {slp:?}"),
    Err(s) => println!("ERROR IN PARSING: \n{s}"),
   } 

    Ok(())
}

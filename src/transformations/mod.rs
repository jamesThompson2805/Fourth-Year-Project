#![allow(unused_imports)]
mod matrix;
mod basis_element;

pub use basis_element::apply_eij_on_metavar;
pub use basis_element::apply_eij_on_program;
pub use basis_element::apply_eijs_on_program;

pub use matrix::apply_eij_poly_on_program;

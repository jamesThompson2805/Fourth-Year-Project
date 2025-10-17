use num_complex::Complex64;

pub type MetaVar = Vec<u32>;

pub type LineVar = usize;

#[derive(Debug, Clone, Copy)]
pub enum SLPVar { L(LineVar), F(Complex64) }

#[derive(Debug, Clone, Copy)]
pub enum Operation { Plus, Mult, }

#[derive(Debug)]
pub enum SLPLine {
    Input(MetaVar),
    Compound((SLPVar, Operation, SLPVar)),
}

pub type SLP = Vec<SLPLine>;
use num_complex::Complex32;

pub type MetaVar = Vec<u32>;

pub type LineVar = u32;

#[derive(Debug)]
pub enum SLPVar { C(MetaVar), L(LineVar), F(Complex32),}

#[derive(Debug, Clone)]
pub enum Operation { Plus, Mult, }

#[derive(Debug)]
pub enum SLPLine {
    Input(SLPVar),
    Compound((SLPVar, Operation, SLPVar)),
}

pub type SLP = Vec<SLPLine>;
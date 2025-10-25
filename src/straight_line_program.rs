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
fn stringify_slpvar(var: &SLPVar) -> String {
    use SLPVar::*;
    match var {
        F(c) => c.to_string(),
        L(lno) => format!("L{lno}"),
    }
}
fn stringify_slpline(line: &SLPLine) -> String {
    use SLPLine::*;
    match line {
        Input(m) => format!("C<{}>",m.into_iter().map(|i| i.to_string()).collect::<Vec<String>>().join(",")),
        Compound((s1,Operation::Plus, s2)) => format!("{} + {}", stringify_slpvar(s1), stringify_slpvar(s2)),
        Compound((s1,Operation::Mult, s2)) => format!("{} * {}", stringify_slpvar(s1), stringify_slpvar(s2)),
    }
}
pub fn stringify_slp(p_slp: &SLP) -> String {
    p_slp.iter().enumerate().map(|(i,l)| "L".to_string() + i.to_string().as_str() + ": " + stringify_slpline(l).as_str() )
        .collect::<Vec<String>>().join("\n")
}

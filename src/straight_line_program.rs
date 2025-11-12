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

/*
    sum_slp takes two slp, adds second to first, adds line to first to sum the slp and wipes second
    does not apply any logic to ensuring no duplication, may duplicate metavariable definition
 */
pub fn add_slp(slp1: &mut SLP, slp2: &mut SLP) {
    use SLPLine::*;
    use SLPVar::*;
    use Operation::*;
    let len1 = slp1.len();
    let len2 = slp2.len();
    slp1.append(slp2);
    slp1.push(
        Compound(( L(len1 - 1), Plus, L(len1 + len2 - 1)))
    );
}
pub fn mult_slp(slp1: &mut SLP, slp2: &mut SLP) {
    use SLPLine::*;
    use SLPVar::*;
    use Operation::*;
    let len1 = slp1.len();
    let len2 = slp2.len();
    slp1.append(slp2);
    slp1.push(
        Compound(( L(len1 - 1), Mult, L(len1 + len2 - 1)))
    );
}

/* 
use std::collections::BTreeMap;
type MetaMonomial = BTreeMap<MetaVar, u64>;
type MetaPolynomial = BTreeMap<MetaMonomial, Complex64>;

fn create_coeff(c: &Complex64) -> MetaPolynomial {
    BTreeMap::from([( BTreeMap::from([(Vec::new(),1)]), *c )])
}

fn evaluate_plus( s1: &SLPVar, s2: &SLPVar, so_far: &Vec<Option<MetaPolynomial>>) -> Option<MetaPolynomial> {
    use SLPVar::*;
    let mp_1 = match s1 {
        L(n) => so_far[*n].clone(),
        F(c) => create_coeff(c),
    };
    let mut mp_2 = match s2 {
        L(n) => so_far[*n].clone(),
        F(c) => create_coeff(c),
    };
    let mut mp: MetaPolynomial = BTreeMap::new();
    for (k,v) in mp_1 {
        if mp_2.contains_key(&k) {
            let v_new = v + mp_2[&k];
            mp_2.remove(&k);
            mp.insert(k, v_new);
        } else {
            mp.insert(k,v);
        }
    }
    for (k,v) in mp_2 {
        mp.insert(k,v);
    }
    unimplemented!()
}

pub fn evaluate_slp(slp: &SLP) -> String {
    use SLPLine::*;
    use Operation::*;
    let mut mp_v: Vec<Option<MetaPolynomial>> = Vec::new();
    for line in slp {
        let mp = match line {
                Input(m) => {
                    let mm:MetaMonomial = BTreeMap::from([(m.clone(),1)]);
                    Some(BTreeMap::from( [(mm, Complex64::ONE)] ))
                },
                Compound((s1,Plus, s2)) => ,
                Compound((s1,Mult, s2)) => ,

        };
        mp_v.push(mp);
        
    }
    unimplemented!()
}
*/
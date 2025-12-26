use std::collections::HashMap;

use num_complex::Complex64;

type Metamonomial = HashMap<Vec<u32>, u32>;
type Metapolynomial = HashMap<Metamonomial, Complex64>;
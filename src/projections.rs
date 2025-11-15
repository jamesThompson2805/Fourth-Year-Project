

fn find_casimir(p: usize, k: usize) -> Vec<Vec<(i,j)>> {
    let mut q = Vec::new();
    q.push(Vec::new());
    let mut cycles = Vec::new();
    while !q.is_empty() {
        let str = q.remove(0);
        for i in 0..k {
            let mut str_clone = str.clone();
            str_clone.push(i);
            if str_clone.len() == p {
                cycles.push(str_clone);
            } else {
                q.push(str_clone);
            }
        }
    }
    cycles.iter_mut().for_each(|&mut str| {str.push( str[0] );});
    cycles.map(|str| {})

}
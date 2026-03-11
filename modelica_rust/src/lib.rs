use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_void};
use std::fs::File;
use std::io::{BufRead, BufReader};

extern "C" {
    fn ModelicaMessage(string: *const c_char);
}

fn modelica_message(msg: &str) {
    let c_str = CString::new(msg).unwrap();
    unsafe {
        ModelicaMessage(c_str.as_ptr());
    }
}

#[derive(Debug, Default, Clone)]
pub struct Taylor {
    pub order: i32,
    pub nrow: i32,
    pub ncol: i32,
    pub nq: i32,
    pub nqn: i32,
    pub structure: i32,
    pub m0: Vec<f64>,
    pub m1: Vec<f64>,
    pub mn: Vec<f64>,
}

#[derive(Debug, Default, Clone)]
pub struct Node {
    pub id: i32,
    pub r_frame: String,
    pub orig: Taylor,
    pub phi: Taylor,
    pub psi: Taylor,
    pub ap: Taylor,
}

#[derive(Debug, Default)]
pub struct SidData {
    pub num_nodes: i32,
    pub num_modes: i32,
    pub use_geo_stiffness: i32,
    pub mass: f64,
    pub nodes: Vec<Node>,
    pub md_cm: Taylor,
    pub j: Taylor,
    pub ct: Taylor,
    pub cr: Taylor,
    pub me: Taylor,
    pub gr: Taylor,
    pub ge: Taylor,
    pub oe: Taylor,
    pub ksigma: Taylor,
    pub ke: Taylor,
    pub de: Taylor,
}

impl Taylor {
    fn new(nrow: i32, ncol: i32, nq: i32, nqn: i32, structure: i32) -> Self {
        Self {
            order: 0,
            nrow,
            ncol,
            nq,
            nqn,
            structure,
            m0: vec![0.0; (nrow * ncol) as usize],
            m1: vec![0.0; (nrow * nq * ncol) as usize],
            mn: if nqn > 0 { vec![0.0; (nrow * nqn * ncol) as usize] } else { Vec::new() },
        }
    }
}

fn get_next_integer(line: &str) -> Option<i32> {
    let mut current_num = String::new();
    let mut found = false;
    for c in line.chars() {
        if c.is_ascii_digit() {
            current_num.push(c);
            found = true;
        } else if found {
            break;
        }
    }
    current_num.parse().ok()
}

// C version of getNextInteger was a bit more complex, it took a startToken.
// Let's implement a version that finds all integers in a line.
fn find_integers(line: &str) -> Vec<i32> {
    let mut results = Vec::new();
    let mut current = String::new();
    for c in line.chars() {
        if c.is_ascii_digit() || c == '-' {
            current.push(c);
        } else {
            if !current.is_empty() {
                if let Ok(val) = current.parse::<i32>() {
                    results.push(val);
                }
                current.clear();
            }
        }
    }
    if !current.is_empty() {
        if let Ok(val) = current.parse::<i32>() {
            results.push(val);
        }
    }
    results
}

fn parse_real(s: &str) -> Option<f64> {
    // SID files use D for exponent
    let s = s.replace('D', "E");
    s.trim().parse().ok()
}

fn find_reals(line: &str) -> Vec<f64> {
    let mut results = Vec::new();
    // This is tricky because reals can have signs and exponents.
    // The C code was very specific.
    // Let's try to split by whitespace and then parse.
    for part in line.split_whitespace() {
        if let Some(val) = parse_real(part) {
            results.push(val);
        }
    }
    results
}

fn parse_taylor<R: BufRead>(reader: &mut R, first_line: &str) -> Taylor {
    let mut order = 0;
    let mut nrow = 0;
    let mut ncol = 0;
    let mut nq = 0;
    let mut nqn = 0;
    let mut structure = 0;
    
    // We might have already read some parts if we follow C logic, but C's parseTaylor starts with fgets.
    // The first line passed here is the one that triggered parseTaylor (e.g. "mdCM")
    
    let mut taylor = Taylor::default();
    let mut initialized = false;

    let mut current_line = first_line.to_string();
    loop {
        if current_line.contains("order ") {
            if let Some(v) = get_next_integer(&current_line) { order = v; }
        } else if current_line.contains("nrow ") {
            if let Some(v) = get_next_integer(&current_line) { nrow = v; }
        } else if current_line.contains("ncol ") {
            if let Some(v) = get_next_integer(&current_line) { ncol = v; }
        } else if current_line.contains("nqn ") {
            if let Some(v) = get_next_integer(&current_line) { nqn = v; }
        } else if current_line.contains("nq ") {
            if let Some(v) = get_next_integer(&current_line) { nq = v; }
        } else if current_line.contains("structure ") {
            if let Some(v) = get_next_integer(&current_line) { structure = v; }
        } else if current_line.contains("m0") {
            if !initialized {
                taylor = Taylor::new(nrow, ncol, nq, nqn, structure);
                taylor.order = order;
                initialized = true;
            }
            let ints = find_integers(&current_line);
            let reals = find_reals(&current_line);
            if ints.len() >= 3 && !reals.is_empty() {
                // m0(i, j) = val. C uses 1-based indexing in file? 
                // C code: idx1 = getNextInteger(line, tokenIdx + 1, &idx1); // idx1
                //         idx2 = getNextInteger(line, tokenIdx + 1, &idx2); // idx2
                //         i = matrix2Index(idx1, idx2, nrow, ncol);
                let r = ints[1];
                let c = ints[2];
                let val = reals[0];
                let idx = ((r - 1) * ncol + (c - 1)) as usize;
                if idx < taylor.m0.len() {
                    taylor.m0[idx] = val;
                }
            }
        } else if current_line.contains("m1") {
            if !initialized {
                taylor = Taylor::new(nrow, ncol, nq, nqn, structure);
                taylor.order = order;
                initialized = true;
            }
            let ints = find_integers(&current_line);
            let reals = find_reals(&current_line);
            if ints.len() >= 4 && !reals.is_empty() {
                let r = ints[1];
                let q = ints[2];
                let c = ints[3];
                let val = reals[0];
                // matrix3Index(int r, int q, int c, int numR, int numQ, int numC)
                // return (q - 1) * numC * numR + (r - 1)*numC + c - 1;
                let idx = ((q - 1) * ncol * nrow + (r - 1) * ncol + (c - 1)) as usize;
                if idx < taylor.m1.len() {
                    taylor.m1[idx] = val;
                }
            }
        } else if current_line.contains("end ") {
            if !initialized {
                taylor = Taylor::new(nrow, ncol, nq, nqn, structure);
                taylor.order = order;
            }
            return taylor;
        }

        current_line.clear();
        if reader.read_line(&mut current_line).unwrap() == 0 { break; }
    }
    taylor
}

#[no_mangle]
pub unsafe extern "C" fn SIDFileConstructor_C(file_name_ptr: *const c_char) -> *mut c_void {
    let file_name = match CStr::from_ptr(file_name_ptr).to_str() {
        Ok(s) => s,
        Err(_) => {
            return std::ptr::null_mut();
        }
    };
    let file = match File::open(file_name) {
        Ok(f) => f,
        Err(_) => {
            return std::ptr::null_mut();
        }
    };

    let mut sid = Box::new(SidData::default());
    let mut reader = BufReader::new(file);
    let mut line = String::new();

    let mut initialized = false;
    while reader.read_line(&mut line).unwrap() > 0 {
        if !initialized {
            let ints = find_integers(&line);
            if ints.len() >= 2 {
                sid.num_nodes = ints[0];
                sid.num_modes = ints[1];
                sid.nodes = vec![Node::default(); sid.num_nodes as usize];
                initialized = true;
            }
            line.clear();
            continue;
        }

        if line.contains("refmod") {
            loop {
                line.clear();
                if reader.read_line(&mut line).unwrap() == 0 { break; }
                if line.contains(" mass ") {
                    let reals = find_reals(&line);
                    if !reals.is_empty() {
                        sid.mass = reals[0];
                    }
                } else if line.contains("end refmod") {
                    break;
                }
            }
        } else if line.contains("new node") {
            let node_id = get_next_integer(&line).unwrap_or(0);
            if let Some(node) = sid.nodes.iter_mut().find(|n| n.id == 0) {
                node.id = node_id;
                loop {
                    line.clear();
                    if reader.read_line(&mut line).unwrap() == 0 { break; }
                    if line.contains("origin") {
                        node.orig = parse_taylor(&mut reader, &line);
                    } else if line.contains("phi") {
                        node.phi = parse_taylor(&mut reader, &line);
                    } else if line.contains("psi") {
                        node.psi = parse_taylor(&mut reader, &line);
                    } else if line.contains("AP") {
                        node.ap = parse_taylor(&mut reader, &line);
                    } else if line.contains("end node") {
                        break;
                    }
                }
            }
        } else if line.contains("mdCM") {
            sid.md_cm = parse_taylor(&mut reader, &line);
        } else if line.contains("J") {
            sid.j = parse_taylor(&mut reader, &line);
        } else if line.contains("Ct") {
            sid.ct = parse_taylor(&mut reader, &line);
        } else if line.contains("Cr") {
            sid.cr = parse_taylor(&mut reader, &line);
        } else if line.contains("Me") {
            sid.me = parse_taylor(&mut reader, &line);
        } else if line.contains("Gr") {
            sid.gr = parse_taylor(&mut reader, &line);
        } else if line.contains("Ge") {
            sid.ge = parse_taylor(&mut reader, &line);
        } else if line.contains("Oe") {
            sid.oe = parse_taylor(&mut reader, &line);
        } else if line.contains("ksigma") {
            sid.ksigma = parse_taylor(&mut reader, &line);
        } else if line.contains("Ke") {
            sid.ke = parse_taylor(&mut reader, &line);
        } else if line.contains("De") {
            sid.de = parse_taylor(&mut reader, &line);
        }
        line.clear();
    }

    Box::into_raw(sid) as *mut c_void
}


#[no_mangle]
pub unsafe extern "C" fn SIDFileDestructor_C(p_sid: *mut c_void) {
    if !p_sid.is_null() {
        drop(Box::from_raw(p_sid as *mut SidData));
    }
}

#[no_mangle]
pub unsafe extern "C" fn getMass(p_sid: *mut c_void) -> f64 {
    if p_sid.is_null() { return 0.0; }
    let sid = &*(p_sid as *mut SidData);
    sid.mass
}

fn get_taylor_by_name<'a>(sid: &'a SidData, name: &str) -> Option<&'a Taylor> {
    match name {
        "mdCM" => Some(&sid.md_cm),
        "J" => Some(&sid.j),
        "Ct" => Some(&sid.ct),
        "Cr" => Some(&sid.cr),
        "Me" => Some(&sid.me),
        "Gr" => Some(&sid.gr),
        "Ge" => Some(&sid.ge),
        "Oe" => Some(&sid.oe),
        "ksigma" => Some(&sid.ksigma),
        "Ke" => Some(&sid.ke),
        "De" => Some(&sid.de),
        _ => None,
    }
}

#[no_mangle]
pub unsafe extern "C" fn getM0(p_sid: *mut c_void, taylor_name_ptr: *const c_char, m0: *mut f64, nr: i32, nc: i32) {
    if p_sid.is_null() { return; }
    let sid = &*(p_sid as *mut SidData);
    let taylor_name = match CStr::from_ptr(taylor_name_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return,
    };
    if let Some(t) = get_taylor_by_name(sid, taylor_name) {
        if nr == t.nrow && nc == t.ncol {
            std::ptr::copy_nonoverlapping(t.m0.as_ptr(), m0, (nr * nc) as usize);
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn getM1(p_sid: *mut c_void, taylor_name_ptr: *const c_char, m1: *mut f64, nr: i32, nq: i32, nc: i32) {
    if p_sid.is_null() { return; }
    let sid = &*(p_sid as *mut SidData);
    let taylor_name = match CStr::from_ptr(taylor_name_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return,
    };
    if let Some(t) = get_taylor_by_name(sid, taylor_name) {
        if nr == t.nrow && nq == t.nq && nc == t.ncol {
            std::ptr::copy_nonoverlapping(t.m1.as_ptr(), m1, (nr * nq * nc) as usize);
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn getM0Node(p_sid: *mut c_void, taylor_name_ptr: *const c_char, node_idx: i32, m0: *mut f64, nr: i32, nc: i32) {
    if p_sid.is_null() { return; }
    let sid = &*(p_sid as *mut SidData);
    let taylor_name = match CStr::from_ptr(taylor_name_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return,
    };
    
    if node_idx > 0 && node_idx <= sid.num_nodes {
        let node = &sid.nodes[(node_idx - 1) as usize];
        let t = match taylor_name {
            "origin" => &node.orig,
            "phi" => &node.phi,
            "psi" => &node.psi,
            "AP" => &node.ap,
            _ => &node.orig,
        };
        if nr == t.nrow && nc == t.ncol {
            std::ptr::copy_nonoverlapping(t.m0.as_ptr(), m0, (nr * nc) as usize);
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn getM1Node(p_sid: *mut c_void, taylor_name_ptr: *const c_char, node_idx: i32, m1: *mut f64, nr: i32, nq: i32, nc: i32) {
    if p_sid.is_null() {
        return;
    }
    let sid = &*(p_sid as *mut SidData);
    let taylor_name = match CStr::from_ptr(taylor_name_ptr).to_str() {
        Ok(s) => s,
        Err(_) => {
            return;
        }
    };
    
    if node_idx > 0 && node_idx <= sid.num_nodes {
        let node = &sid.nodes[(node_idx - 1) as usize];
        let t = match taylor_name {
            "origin" => &node.orig,
            "phi" => &node.phi,
            "psi" => &node.psi,
            "AP" => &node.ap,
            _ => {
                &node.orig
            }
        };
        
        let expected_size = (nr as i64) * (nq as i64) * (nc as i64);
        if nr == t.nrow && nq == t.nq && nc == t.ncol {
            if (t.m1.len() as i64) >= expected_size {
                std::ptr::copy_nonoverlapping(t.m1.as_ptr(), m1, expected_size as usize);
            }
        }
    }
}


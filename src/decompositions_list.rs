#![allow(dead_code)]

use num::*;
use quizx::{scalar::{ScalarN, Scalar, FromPhase, Sqrt2}, hash_graph::{GraphLike, VType, EType}};
use quizx::vec_graph::Graph;
use array_tool::vec::Intersect;
use std::cmp::{Eq, Ord, Reverse};
use std::collections::{BinaryHeap, HashMap};
use std::hash::Hash;


/// Pick the first <= 6 T gates from the given graph
pub fn first_ts(g: &Graph) -> Vec<usize> {
    let mut t = vec![];

    for v in g.vertices() {
        if *g.phase(v).denom() == 4 { t.push(v); }
        if t.len() == 6 { break; }
    }

    t
}

//return the 6 most connected T gates
pub fn most_connected_ts(g :&Graph) -> Vec<usize>{

    let ts:Vec<usize> = g.vertices().filter(|&v| *g.phase(v).denom() == 4 ).collect();
    let mut best_ts = vec![];
    for _ in 0..6{
        let best = ts.iter().filter(|&v| !best_ts.contains(v)).max_by_key(|&v| g.neighbors(*v).len()).expect("this shouldn't be printed");
        best_ts.push(*best);
    }

    best_ts
}


fn replace_t0(g:  &Graph, verts: &[usize]) -> Graph {
    // println!("replace_t0");
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);
    let w = g.add_vertex(VType::Z);
    g.add_edge_with_type(verts[0], w, EType::H);
    g.add_to_phase(verts[0], Rational::new(-1,4));
    g
}

fn replace_t1(g: &Graph, verts: &[usize]) -> Graph {
    // println!("replace_t1");
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![1,0,1,0]);
    let w = g.add_vertex_with_phase(VType::Z, Rational::one());
    g.add_edge_with_type(verts[0], w, EType::H);
    g.add_to_phase(verts[0], Rational::new(-1,4));
    g
}

pub fn trivial_replace(g: &Graph, verts: &[usize])-> Vec<Graph>{
    vec![replace_t0(g,verts),replace_t1(g,verts)]
}



fn replace_magic5_0(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(1, vec![1, 0, 0, 0]);    
    for &v in verts {
        g.add_to_phase(v, Rational::new(-1,4));
        g.add_edge_smart(v, verts[0], EType::N);
    }
    g.add_to_phase(verts[0], Rational::new(-3,4));
    g
}

fn replace_magic5_1(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();  
    *g.scalar_mut() *= ScalarN::Exact(1, vec![-1, 0, 1, 0]);
    let p = g.add_vertex(VType::Z);
    for &v in verts {
        g.add_to_phase(v, Rational::new(-1,4));
        g.add_edge_with_type(v, p, EType::H);
    }
    let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1,4));
    g.add_edge_with_type(w, p, EType::H);
    g
}

fn replace_magic5_2(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(9, vec![0, -1, 0, 0]);
    let p = g.add_vertex(VType::Z);
    let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1,4));
    g.add_edge_with_type(p, w, EType::H);
    for i in 0..verts.len() {
        g.add_to_phase(verts[i], Rational::new(-1,4));
        g.add_edge_with_type(verts[i], p, EType::H);
        g.add_edge_with_type(verts[i], w, EType::H);
        for j in i+1..verts.len() {
            g.add_edge_smart(verts[i], verts[j], EType::H);
        }
    }
    g
}

pub fn replace_magic5(g: &Graph, verts: &[usize])-> Vec<Graph>{
    vec![replace_magic5_0(g,verts),replace_magic5_1(g,verts),replace_magic5_2(g,verts)]
}


/// Returns a best occurrence of a cat state
/// The fist vertex in the result is the Clifford spider
pub fn cat_ts(g: &Graph) -> Vec<usize> {
    // the graph g is supposed to be completely simplified

    // g.vertices()
    //     .filter(|&v| g.phase(v).is_integer() && g.degree(v) <= 6)
    //     .min_by_key(|v| [4, 6, 5, 3].into_iter()
    //         .position(|r| r == g.degree(r)))
    //     .map(|v| {
    //         let mut neighbors = g.neighbor_vec(v);
    //         neighbors.insert(0, v);
    //         neighbors
    //     })
    //     .unwrap_or(Vec::new())


    let prefered_order = [4,6,5,3];

    let mut res = vec![];
    let mut index = None;
    for v in g.vertices() {
        if g.phase(v).denom() == &1{
            let mut neigh = g.neighbor_vec(v);
            if neigh.len() <= 6 {
                match prefered_order.iter().position(|&r| r == neigh.len()) {
                    Some(this_ind) => match index {
                        Some(ind) if this_ind < ind => {res = vec![v]; res.append(&mut neigh); index = Some(this_ind);},
                        None => {res = vec![v]; res.append(&mut neigh); index = Some(this_ind);},
                        _ => (),
                    },
                    _ => (),
                }
                if index == Some(0){break;}
            }
        }
    }
    res
}

//return a cat state with the best arity which is also the most connected of it's type 
pub fn cat_ts_connected_first(g: &Graph) -> Vec<usize> {
    // the graph g is supposed to be completely simplified
    let prefered_order = [4,6,5,3];
    for k in prefered_order{
        if count_catk(g,k)!=0{
            let catsk:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer() && g.degree(v) == k).collect();
            let mut best_cat = 0;
            let mut best_arity = 0;
            for c in catsk{
                let n_cat_arity = cat_arity_for_cat(&g, c);
                if n_cat_arity > best_arity {
                    best_cat = c;
                    best_arity = n_cat_arity;
                }
            }
            let mut best_full_cat = vec![best_cat];
            best_full_cat.append(&mut g.neighbor_vec(best_cat));
            return best_full_cat;
        }   
    }
    vec![]
}



pub fn replace_cat_decomp(g: &Graph, verts: &[usize]) -> Vec<Graph> {
    // verts[0] is a 0- or pi-spider, linked to all and only to vs in verts[1..] which are T-spiders
    let mut g = g.clone(); // that is annoying ... 
    let mut verts = Vec::from(verts);
    if g.phase(verts[0]).numer() == &1 {
        g.set_phase(verts[0], Rational::new(0,1));
        let mut neigh = g.neighbor_vec(verts[1]);
        neigh.retain(|&x| x != verts[0]);
        for &v in &neigh{
            g.add_to_phase(v, Rational::new(1,1));
        }
        let tmp = g.phase(verts[1]);
        *g.scalar_mut() *= ScalarN::from_phase(tmp);
        g.set_phase(verts[1], g.phase(verts[1])*Rational::new(-1,1));
    }
    if [3,5].contains(&verts[1..].len()) {
        let w = g.add_vertex(VType::Z);
        let v = g.add_vertex(VType::Z);
        g.add_edge_with_type(v, w, EType::H);
        g.add_edge_with_type(v, verts[0], EType::H);
        verts.push(v);     
    }
    if verts[1..].len() == 6 {
        return vec![replace_cat6_0(&g,&verts),replace_cat6_1(&g,&verts),replace_cat6_2(&g,&verts)]

    }else if verts[1..].len() == 4 {
       return vec![replace_cat4_0(&g, &verts),replace_cat4_1(&g, &verts)]
    }else { println!("this shouldn't be printed"); };

    vec![g]
}

fn replace_cat4_0(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(0, vec![0, 0, 1, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
    }
    g
}

fn replace_cat4_1(g: &Graph, verts: &[usize]) -> Graph {
    // same as replace_cat6_0, only with a different scalar
    let mut g = g.clone();  
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, -1, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
        g.set_edge_type(v, verts[0], EType::N);
    }
    g.set_phase(verts[0], Rational::new(-1,2));
    g
}


fn replace_cat6_0(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, 0, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
        g.set_edge_type(v, verts[0], EType::N);
    }
    g.set_phase(verts[0], Rational::new(-1,2));
    g
}

fn replace_cat6_1(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();  
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![-1, 0, 1, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
    }
    g
}

fn replace_cat6_2(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(7, vec![0, -1, 0, 0]);
    for i in 1..verts.len() {
        g.add_to_phase(verts[i], Rational::new(-1,4));
        for j in i+1..verts.len() {
            g.add_edge_smart(verts[i], verts[j], EType::H);
        }
    }
    g
}

pub fn find_2cat3(g :&Graph) -> Vec<usize>{
    //TODO make it work with pi phase on cat g.phase(v).is_integer()
    let cats3:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_zero() && g.degree(v) == 3).collect();
    
    if cats3.len() < 2 {
        return vec![];
    }
    
    for i in 0..cats3.len()-1 {
        for &c2 in &cats3[i+1 ..] {

           let n1 = g.neighbor_vec(cats3[i]);
           let n2 = g.neighbor_vec(c2);

            if n1.intersect(n2.clone()).len() == 0 {

                let mut tab = vec![cats3[i],c2];
                tab.extend(n1);
                tab.extend(n2);
                return tab
            }
        }
    }
    vec![]
}

pub fn replace_2cat3(g :&Graph,verts: &[usize])-> Vec<Graph> {
    vec![replace_2cat3_0(&g,&verts),replace_2cat3_1(&g,&verts),replace_2cat3_2(&g,&verts)]
}


fn replace_2cat3_0(g :&Graph,verts: &[usize]) -> Graph{

    let mut g = g.clone();
    
    //(1+i)/2
    *g.scalar_mut() *= ScalarN::from_phase(Rational::new(1,4))* Scalar::sqrt2_pow(-1);

    g.add_to_phase(verts[2], Rational::new(1,4));
    g.add_to_phase(verts[3], Rational::new(-1,4));
    g.add_to_phase(verts[4], Rational::new(1,4));
    g.add_to_phase(verts[5], Rational::new(-1,4));
    g.add_to_phase(verts[6], Rational::new(-1,4));
    g.add_to_phase(verts[7], Rational::new(-1,4));
    //quizx::simplify::full_simp(&mut g);
    g
}


fn replace_2cat3_1(g :&Graph,verts: &[usize]) -> Graph{

    let mut g = g.clone();
    
    //(1 - 1j)/2**0.5
    *g.scalar_mut() *= ScalarN::from_phase(Rational::new(-1,4));

    g.add_to_phase(verts[2], Rational::new(-1,4));
    g.add_to_phase(verts[3], Rational::new(1,4));
    g.add_to_phase(verts[4], Rational::new(-1,4));
    g.add_to_phase(verts[5], Rational::new(3,4));
    g.add_to_phase(verts[6], Rational::new(-1,4));
    g.add_to_phase(verts[7], Rational::new(-1,4));
    g.add_edge_smart(verts[6], verts[7], EType::H);
    //quizx::simplify::full_simp(&mut g);
    g
}

fn replace_2cat3_2(g :&Graph,verts: &[usize]) -> Graph{

    let mut g = g.clone();
    
    //(0.5 + i)
    *g.scalar_mut() *= ScalarN::from_phase(Rational::new(1,2));
    *g.scalar_mut() *= ScalarN::sqrt2_pow(-6);

    for i in 2..8{
        g.add_to_phase(verts[i], Rational::new(-1,4));
        
        let nv = g.add_vertex(VType::X);
        g.add_edge(verts[i], nv);
        if i==2 || i==4 {
            g.add_to_phase(nv,  Rational::new(1,1));
        }
    }
    g.remove_vertex(verts[0]);
    g.remove_vertex(verts[1]);
    //quizx::simplify::full_simp(&mut g);
    g
}


pub fn cat_distribution(g: &Graph) -> Vec<usize>{


    let mut cats:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer()).map(|x| g.neighbors(x).len()).collect();

    cats.sort();

    cats
}

pub fn count_2cat3(g :&Graph,intersection : usize) -> usize{

    let mut counter = 0;

    let cats3:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer() && g.degree(v) == 3).collect();
    for i in 0..cats3.len()-1 {
        for &c2 in &cats3[i+1 ..] {

           let n1 = g.neighbor_vec(cats3[i]);
           let n2 = g.neighbor_vec(c2);

            if n1.intersect(n2.clone()).len() == intersection {

                counter += 1;
            }
        }
    }
    counter
}

pub fn count_3_connected_cat3(g :&Graph,intersection : usize) -> usize{

    let mut counter = 0;

    let cats3:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer() && g.degree(v) == 3).collect();
    for i in 0..cats3.len()-2 {
        for j in i..cats3.len()-1 {
            for k in j..cats3.len() {

                let n1 = g.neighbor_vec(cats3[i]);
                let n2 = g.neighbor_vec(cats3[j]);
                let n3 = g.neighbor_vec(cats3[k]);

                if n1.intersect(n2.clone()).intersect(n3.clone()).len() == intersection {

                    counter += 1;
                }
            }
        }
    }
    counter
}

pub fn most_connected_t_incat3(g :&Graph) -> (usize,usize){

    let cats3:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer() && g.degree(v) == 3).collect();
    if cats3.len() == 0 {return (0,0)}
    let mut ts = vec![];

    for c in cats3 {
        ts.append(&mut g.neighbor_vec(c));
    }

    let (count,x) = most_frequent(&ts, 1);

    //println!("best t gate is in {} cat 3",count);
    (*x,count)
}

pub fn most_connected_t_incatk(g :&Graph,k:usize) -> (usize,usize){

    let catsk:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer() && g.degree(v) == k).collect();
    if catsk.len() == 0 {return (0,0)}
    let mut ts = vec![];

    for c in catsk {
        ts.append(&mut g.neighbor_vec(c));
    }

    let (count,x) = most_frequent(&ts, 1);

    //println!("best t gate is in {} cat 3",count);
    (*x,count)
}












// from https://stackoverflow.com/questions/64262297/rust-how-to-find-n-th-most-frequent-element-in-a-collection
fn most_frequent<T>(array: &[T], k: usize) -> (usize, &T) where T: Hash + Eq + Ord,{
    let mut map = HashMap::with_capacity(array.len());
    for x in array {
        *map.entry(x).or_default() += 1;
    }

    let mut heap = BinaryHeap::with_capacity(k + 1);
    for (x, count) in map.into_iter() {
        if heap.len() < k {
            heap.push(Reverse((count, x)));
        } else {
            let &Reverse((min, _)) = heap.peek().unwrap();
            if min < count {
                heap.pop();
                heap.push(Reverse((count, x)));
            }
        }
    }
    let result:Vec<(usize, &T)> = heap.into_sorted_vec().into_iter().map(|r| r.0).collect();
    result[k - 1]
}


pub fn count_catk(g :&Graph,k : usize) -> usize{

    let catsk:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer() && g.degree(v) == k).collect();
    catsk.len()
}

//return the number of cat-k center vertices neiboring a vertex 
fn cat_k_arity(g :&Graph,v:usize,k : usize) -> usize{

    let catsk:Vec<usize> = g.neighbors(v).filter(|&v| g.phase(v).is_integer() && g.degree(v) == k).collect();
    catsk.len()
}

//return the number of cat connected to a cat centered at k
fn cat_arity_for_cat(g :&Graph,v:usize) -> usize {
    let mut cats = vec![];

//    .filter(|&v| g.phase(v).is_integer() && g.degree(v) == k).collect();

    for t in g.neighbors(v) {
        cats.append(&mut g.neighbors(t).filter(|&x| g.phase(x).is_integer() && v!=x).collect());   
    }

    cats.sort();
    cats.dedup();
    cats.len()
}


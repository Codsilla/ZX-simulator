
#[allow(unreachable_code)]

use num::{Rational, One};
use quizx::{vec_graph::Graph, hash_graph::{GraphLike, VType}, scalar::*};


pub fn cutter(g :& Graph,cuts :& Vec<usize>) -> Vec<Graph> {

    recursive_cut(&mut cuts.clone(), vec![g.clone()] ).1
}

fn recursive_cut(cuts :&mut Vec<usize>, stack : Vec<Graph>) -> (&mut Vec<usize>, Vec<Graph>){

    if cuts.len() == 0 {
        return (cuts,stack)
    }

    let c = cuts.pop().expect("Cant pop an empty array of cuts");
    let mut nstack = Vec::<Graph>::new();

    for mut g in stack {

        let mut g2 = g.clone();

        //normalization for implicitly expressing |0>,|1> as xspider (1 /sqrt(2))
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);
        *g2.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);


        
        for n in 0..g2.neighbors(c).len() {
            let neigh = g2.neighbor_at(c, n);
            g2.add_to_phase(neigh,Rational::new(1,1));

            //normalization for implicitly using the pi-copy rule (1 /sqrt(2)) for n-1 neighbors
            if n!=0 {
                *g.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);
                *g2.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);
            }
        }



        //normalization for the angle in the spider removed
        let phase = g2.phase(c);
        g2.scalar_mut().mul_phase(phase);

        g.remove_vertex(c);
        g2.remove_vertex(c);

        nstack.push(g);
        nstack.push(g2);

    }

    recursive_cut(cuts,nstack)
}


pub fn to_subcomponent(g :& Graph) -> Vec<Graph>{


    let comp = g.component_vertices();
    let mut subzx = vec![g.clone();comp.len()];
    

    for i in 0..comp.len(){    
        for v in g.vertices() {

            if !comp[i].contains(&v) {
                subzx[i].remove_vertex(v);
            }
        }

        //remove phase except for the first one
        if i!=0 { 
            subzx[i].scalar_mut().set_one();
            
        }
    }



    subzx
}


pub fn is_cut_worth(g :& Graph,cuts :& Vec<usize>, alpha : f64) -> bool {
    
    if cuts.len()<1 || g.tcount() < 24 { return false}
    let mut g = g.clone();
    let t_before = g.tcount();



    //do the cut
    for &v in cuts {
        g.remove_vertex(v);
    }

    quizx::simplify::full_simp(&mut g);

    let comp = g.component_vertices();

    if comp.len() == 0 { return (cuts.len() as f64) < alpha * (t_before as f64); }

    let mut tcounts : Vec<isize> = Vec::new();
    //count number of T in each component by removing them one by one
    for mut c in comp {
        c.sort();
        c.dedup();
        let t = g.tcount();
        

        for i in 0..c.len() {
            let v = c[i];
            g.remove_vertex(v);
        }

        tcounts.push(t as isize - g.tcount() as isize);
    }
    let t_after = *tcounts.iter().max().expect("something went wrong");

    cuts.len() as f64 + alpha*(t_after as f64) + 1.0 < alpha * (t_before as f64)
}


pub fn direct_cut(g :& Graph,cuts :& Vec<usize>) -> Vec<Graph> {

    let mut tab = vec![g.clone()];

    for c in cuts {
        let mut newtab = Vec::<Graph>::new();

        for mut g in tab {

            *g.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);
            let nv = g.add_vertex(VType::X);
            g.add_edge(*c, nv);
            let mut g2 = g.clone(); 
            g2.set_phase(nv, Rational::new(1,1));

            newtab.push(g);
            newtab.push(g2);

        }

        tab = newtab;
    }

    tab
}

pub fn cut_index(g :& Graph,cuts :& Vec<usize>, index : usize) -> Graph {

    let mut ng = g.clone();
    let mut i = 0;

    for c in cuts {

        //normalization for basis state
        *ng.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);

        let nv = ng.add_vertex(VType::X);
        ng.add_edge(*c, nv);


        if index >> i & 1 == 1 {
            ng.set_phase(nv, Rational::new(1,1));
        }
        i += 1;
    }
    quizx::simplify::full_simp(&mut ng);
    ng
}

//compute the "effective" alpha relate to a cut
pub fn cut_alpha(g :& Graph,cuts :& Vec<usize>) -> f64 {

    let t_remove = g.tcount() - t_after_cut(&g, &cuts); 

    (cuts.len() as f64) / (t_remove as f64)

}


fn t_after_cut(g :& Graph,cuts :& Vec<usize>) -> usize{

    let mut  g = g.clone();
    
    for &v in cuts {
        g.remove_vertex(v);
    }

    quizx::simplify::full_simp(&mut g);

    let comp = g.component_vertices();

    if comp.len() == 0 { return 0 }

    let mut tcounts : Vec<isize> = Vec::new();
    //count number of T in each component by removing them one by one
    for mut c in comp {
        c.sort();
        c.dedup();
        let t = g.tcount();
        

        for i in 0..c.len() {
            let v = c[i];
            g.remove_vertex(v);
        }

        tcounts.push(t as isize - g.tcount() as isize);
    }
    let t_after = *tcounts.iter().max().expect("something went wrong");

    t_after as usize
}



pub fn effective_alpha_decomp(g: &Graph,decomp: &Vec<Graph>) -> f64 {

    let nb_term = decomp.len();
    let mut t_counts = vec![];

    for g in decomp {

        let mut d = g.clone();
        let comp = d.component_vertices();
        if comp.len() == 0 { continue; }

        for mut c in comp {
            c.sort();
            c.dedup();
            let t = d.tcount();
            
    
            for i in 0..c.len() {
                let v = c[i];
                d.remove_vertex(v);
            }
    
            t_counts.push(t as isize - d.tcount() as isize);
        }
    }
    let t_after = *t_counts.iter().max().expect("something went wrong") as usize;
    
    effective_alpha_by_tcount(g.tcount(),t_after,nb_term)
}

fn effective_alpha_by_tcount(t_before : usize,t_after : usize, nb_term: usize) -> f64 {

    let t_remove = t_before - t_after;
    return (t_remove as f64).log2()/(nb_term as f64)
}

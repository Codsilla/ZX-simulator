#![allow(non_snake_case)]

use crate::{
    cutter::{
        cut_alpha, cut_index, direct_cut, effective_alpha_decomp, is_cut_worth, to_subcomponent,
    },
    decompositions_list::{
        cat_ts, cat_ts_connected_first, find_2cat3, first_ts, most_connected_t_incat3,
        most_connected_ts, replace_2cat3, replace_cat_decomp, replace_magic5
    },
    hypergraph_builder::kaHyPar_cut_finder,
};
use num::*;
use quizx::vec_graph::Graph;
use quizx::{
    decompose::Decomposer,
    hash_graph::GraphLike,
    scalar::{Scalar, ScalarN},
};
use rayon::prelude::*;

//(mut g : Graph, layer : usize, bar : &ProgressBar,imbalance : f64)
fn recursive_decomp_internal(
    mut g: Graph,
    layer: usize,
    imbalance: f64,
    connected_cat3: bool,
    most_connected_cat: bool,
    most_connected_t: bool,
    catcat3: bool,
) -> Scalar<Vec<isize>> {
    quizx::simplify::full_simp(&mut g);

    if g.tcount() < 30 {
        let mut d = Decomposer::new(&g);
        d.use_cats(true);
        d.with_full_simp();
        d.decomp_all();
        return d.scalar;
    }

    let cutKaHypar = kaHyPar_cut_finder(&g, imbalance);
    let (best_t_cat3, count) = most_connected_t_incat3(&g);
    let cutConnectedCat3 = if count > 0 { vec![best_t_cat3] } else { vec![] };
    let cut = if connected_cat3 && (cut_alpha(&g, &cutKaHypar) > cut_alpha(&g, &cutConnectedCat3)) {
        cutConnectedCat3
    } else {
        cutKaHypar
    };
    let nbterm_after_cut = 1 << cut.len();

    if cut.len() != 0 && is_cut_worth(&g, &cut, 0.24) {
        //println!("cut was worth");
        //add the scalars of all diagrams after the cut

        (0..nbterm_after_cut)
            .into_par_iter()
            .map(|i| {
                if layer == 0 { //&& i%((nbterm_after_cut as f64/100.0).ceil() as usize) == 0 {
                     //print!("({}/{})",i,nbterm_after_cut);
                     //io::stdout().flush().unwrap(); //not cool

                    // bar.inc(1);
                }

                let z = cut_index(&g, &cut, i as usize);

                if z.component_vertices().len() == 0 {
                    return z.scalar().clone();
                }

                let mut local_result = ScalarN::one();

                //compute scalar of subdiagrams and multiply them
                for mut subzx in to_subcomponent(&z).into_iter() {
                    quizx::simplify::full_simp(&mut subzx);
                    let plz = recursive_decomp_internal(
                        subzx,
                        layer + 1,
                        imbalance,
                        connected_cat3,
                        most_connected_cat,
                        most_connected_t,
                        catcat3,
                    ); //recursive_decomp_internal(subzx,layer+1,&bar,imbalance);
                       //println!("{}",plz);

                    local_result *= plz;
                }

                //result = result + local_result;
                local_result
            })
            .reduce(|| ScalarN::zero(), |a, b| a + b)

        //return result
    } else {
        let catcat3s = find_2cat3(&g);
        let (x, count) = most_connected_t_incat3(&g);
        let mut result = ScalarN::zero();
        let best_cat = if !most_connected_cat {
            cat_ts(&g)
        } else {
            cat_ts_connected_first(&g)
        };

        if catcat3 && catcat3s.len() != 0 && best_cat.len() !=5 && (!connected_cat3 || count < 2) {
            for d in replace_2cat3(&g, &catcat3s) {
                result = result
                    + recursive_decomp_internal(
                        d,
                        layer + 1,
                        imbalance,
                        connected_cat3,
                        most_connected_cat,
                        most_connected_t,
                        catcat3,
                    ) //recursive_decomp_internal(d,layer+1,&bar,imbalance)
            }
        } else if (count > 1) && connected_cat3 {
            //println!("got here");
            for d in direct_cut(&g, &vec![x]) {
                result = result
                    + recursive_decomp_internal(
                        d,
                        layer + 1,
                        imbalance,
                        connected_cat3,
                        most_connected_cat,
                        most_connected_t,
                        catcat3,
                    ) //recursive_decomp_internal(d,layer+1,&bar,imbalance)
            }
        } else if best_cat.len() != 0 {
            for d in replace_cat_decomp(&g, &best_cat) {
                result = result
                    + recursive_decomp_internal(
                        d,
                        layer + 1,
                        imbalance,
                        connected_cat3,
                        most_connected_cat,
                        most_connected_t,
                        catcat3,
                    ) //recursive_decomp_internal(d,layer+1,&bar,imbalance)
            }
        } else {
            let first_ts = first_ts(&g);
            let most_connected_ts = most_connected_ts(&g);
            let mut ts2replace = first_ts.clone();

            if most_connected_t {
                let score_first = effective_alpha_decomp(&g, &replace_magic5(&g, &first_ts[..5]));
                let score_most_connected =
                    effective_alpha_decomp(&g, &replace_magic5(&g, &most_connected_ts[..5]));
                if score_first > score_most_connected {
                    ts2replace = most_connected_ts;
                }
            };

            for d in replace_magic5(&g, &ts2replace[..5]) {
                result = result
                    + recursive_decomp_internal(
                        d,
                        layer + 1,
                        imbalance,
                        connected_cat3,
                        most_connected_cat,
                        most_connected_t,
                        catcat3,
                    ) //recursive_decomp_internal(d,layer+1,&bar,imbalance)
            }
        }

        return result;
    }
}

pub fn recursive_decomp(g: &Graph, imbalance: f64) -> Scalar<Vec<isize>> {
    //let imbalance = 0.4;
    let g = g.clone();
    //let cut =  kaHyPar_cut_finder(&g,imbalance);
    //let nbterm_after_cut = 1<<cut.len();
    //println!("intial cut size : {}",cut.len());
    //let bar = ProgressBar::new(nbterm_after_cut);
    let x = recursive_decomp_internal(g, 0, imbalance, true, false, false, false); //recursive_decomp_internal(g, 0, &bar,imbalance);
                                                                                   //bar.finish();
    x
}

pub fn recursive_decomp_connected_cat3(
    g: &Graph,
    imbalance: f64,
    connected_cat3: bool,
) -> Scalar<Vec<isize>> {
    //let imbalance = 0.4;
    let g = g.clone();
    //let cut =  kaHyPar_cut_finder(&g,imbalance);
    //let nbterm_after_cut = 1<<cut.len();
    //println!("intial cut size : {}",cut.len());
    //let bar = ProgressBar::new(nbterm_after_cut);
    let x = recursive_decomp_internal(g, 0, imbalance, connected_cat3, false, false, false); //recursive_decomp_internal(g, 0, &bar,imbalance);
                                                                                             //bar.finish();
    x
}

pub fn decomp_most_connected_cat(g: &Graph, imbalance: f64) -> Scalar<Vec<isize>> {
    let g = g.clone();
    recursive_decomp_internal(g, 0, imbalance, true, true, false, false)
}

pub fn decomp_most_connected_t(g: &Graph, imbalance: f64) -> Scalar<Vec<isize>> {
    let g = g.clone();
    recursive_decomp_internal(g, 0, imbalance, true, true, true, false)
}

pub fn decomp_2cat3(g: &Graph, catcat3: bool) -> Scalar<Vec<isize>> {
    let g = g.clone();
    recursive_decomp_internal(g, 0, 0.5, true, true, true, catcat3)
}

pub fn decomp_2cat3_without_most_connected_cat(g: &Graph, catcat3: bool) -> Scalar<Vec<isize>> {
    let g = g.clone();
    recursive_decomp_internal(g, 0, 0.5, true, true, true, catcat3)
}

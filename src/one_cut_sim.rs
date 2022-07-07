use indicatif::ProgressBar;
use quizx::decompose::Decomposer;
use quizx::hash_graph::GraphLike;
use quizx::vec_graph::Graph;
use quizx::scalar::*;
use crate::cutter::{to_subcomponent, is_cut_worth, direct_cut, cut_index};
use crate::hypergraph_builder::{kaHyPar_cut_finder};




pub fn one_cut_sim(g :&Graph)-> Scalar<Vec<isize>> {

    let cut = kaHyPar_cut_finder(&g,0.3);
    println!("cut size {}",&cut.len());
    
    let nbterm_after_cut = 1<<cut.len();
    let mut result = ScalarN::zero();    
    let bar = ProgressBar::new(nbterm_after_cut);

    //add the scalars of all diagrams after the cut
    for z in direct_cut(&g,&cut){
        
        if z.component_vertices().len() == 0  {
            result = result + z.scalar();
            continue;
        }  

        let mut local_result = ScalarN::one();
        
        //compute scalar of subdiagrams and multiply them
        for mut subzx in to_subcomponent(&z).into_iter(){


            quizx::simplify::full_simp(&mut subzx);

            let mut d = Decomposer::new(&subzx);
            d.use_cats(true);
            d.with_full_simp();
            let d = d.decomp_parallel(3);
            local_result *= d.scalar;

        } 
        
        result = result + local_result;
        bar.inc(1);
    }

    bar.finish();
    result
}


pub fn one_cut_sim_low_memory(g :&Graph)-> Scalar<Vec<isize>>{

    let cut = kaHyPar_cut_finder(&g,0.3);
    println!("cut size {}",&cut.len());
    
    let nbterm_after_cut = 1<<cut.len();
    let mut result = ScalarN::zero();    
    let bar = ProgressBar::new(nbterm_after_cut);

    //add the scalars of all diagrams after the cut
    for i in 0..nbterm_after_cut {
        
        let z = cut_index(&g,&cut,i as usize); 

        if z.component_vertices().len() == 0  {
            result = result + z.scalar();
            continue;
        }  

        let mut local_result = ScalarN::one();
        
        //compute scalar of subdiagrams and multiply them
        for mut subzx in to_subcomponent(&z).into_iter(){


            quizx::simplify::full_simp(&mut subzx);

            let mut d = Decomposer::new(&subzx);
            d.use_cats(true);
            d.with_full_simp();
            let d = d.decomp_parallel(3);
            local_result *= d.scalar;

        } 
        
        result = result + local_result;
        bar.inc(1);
    }

    bar.finish();
    result

}



pub fn cut_then_decompose(g :&Graph)-> Scalar<Vec<isize>> {

    let mut g = g.clone();
    let cut =  if g.tcount() > 24 {kaHyPar_cut_finder(&g,0.4)}else{Vec::new()};

    let nbterm_after_cut = 1<<cut.len();
    let mut result = ScalarN::zero();   

    
    if cut.len() !=0 && is_cut_worth(&g,&cut,0.24) {
        //println!("cut was worth");
        //add the scalars of all diagrams after the cut
        for i in 0..nbterm_after_cut {
        
            let z = cut_index(&g,&cut,i as usize); 

            if z.component_vertices().len() == 0  {
                result = result + z.scalar();
                continue;
            }  
    
            let mut local_result = ScalarN::one();
            
            //compute scalar of subdiagrams and multiply them
            for mut subzx in to_subcomponent(&z).into_iter(){

                quizx::simplify::full_simp(&mut subzx);
                let plz = cut_then_decompose(&subzx);
                //println!("{}",plz);

                local_result *= plz;
                
    
            } 
            
            result = result + local_result;
        }
    } else {

        quizx::simplify::full_simp(&mut g);

        let mut d = Decomposer::new(&g);
        d.use_cats(true);
        d.with_full_simp();
        let d = d.decomp_parallel(3);
        result = result + d.scalar;


    }



    result
}
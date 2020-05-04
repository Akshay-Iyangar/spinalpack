#include "localopt.h"

//compute cut_off point
double find_cut_off(array2<double>* cur_mat_p, int cut_offpnt){
    array<double> all_cur(((*cur_mat_p).high1()+1) * ((*cur_mat_p).high2()+1));
    int all_cur_in = 0;
    int in1, in2;
    for (in1 = 0; in1 <= (*cur_mat_p).high1(); in1++ )
	for(in2= 0; in2 <= (*cur_mat_p).high2(); in2++){
		all_cur[all_cur_in] = (*cur_mat_p)(in1, in2);
		all_cur_in ++;
    	}
    all_cur.sort();
    return(all_cur[cut_offpnt]);
}

//used for (reverse) sorting. useful in local pack. 
int comp_triple(const three_tuple<int,int,double>& n_t1, const three_tuple<int,int,double>& n_t2){
    if (n_t1.third() > n_t2.third()) return -1;
    else if (n_t1.third() < n_t2.third()) return 1;
    else return 0;
}

void optimize_local_score(array2<double>* cur_matrix_p, array2<double>* simmatrix_p, 
			  graph ppi1, graph ppi2, double alpha_cons,
			  int iter_cnt, int cut_offpnt, bool greedy_choice){
    int len1 = ppi1.number_of_nodes();
    int len2 = ppi2.number_of_nodes();
 
    array2<double> prev_matrix(len1, len2);
    int in1,in2;
    for (in1=0; in1<len1; in1++)
	for (in2=0; in2<len2; in2++)
		prev_matrix(in1,in2) = (*cur_matrix_p)(in1, in2);

    //create deg_from_index array for both ppis
    node v1, v2, neigh1, neigh2;
    array<node> in_to_node1(len1);
    array<node> in_to_node2(len2); 
    forall_nodes(v1, ppi1) in_to_node1[index(v1)] = v1;
    forall_nodes(v2, ppi2) in_to_node2[index(v2)] = v2;
 
    cut_offpnt = len1 * len2 - cut_offpnt;
    double cut_off; 
    double g_sim; //g_sim indicates pairwise similarity

    //local match variables
    dictionary< int,bool > l_m1, l_m2; //local match1 and 2 
    three_tuple<int,int,double> n_t;
    list< three_tuple<int, int, double> > tuple_l;

    //optimum variables
    edge e;	
    graph local_bg;
    node temp1, temp2;
    list<node> part1, part2;
    edge_array<double> weight;
    map<node, node> from_ppi_to_part1;
    map<node, node> from_part_to_ppi1;
    map<node, node> from_ppi_to_part2;
    map<node, node> from_part_to_ppi2;
    list<edge> match_list;
    //main loop of iterations for local score optimization
    for (int iter_index=1; iter_index <= iter_cnt; iter_index++){
        cut_off = find_cut_off(&prev_matrix, cut_offpnt);
	//1- apply local set packing 
	forall_nodes(v1, ppi1){
	     in1 = index(v1); 
	     if (in1 == len1 / 2)
  	    	  printf("%d\% Progress...\n", iter_index * 100 / iter_cnt - (50 / iter_cnt));
	     if (in1 == len1 - 1)
  	    	  printf("%d\% Progress...\n", iter_index * 100 / iter_cnt);
	     forall_nodes(v2, ppi2){
		  in2 = index(v2);
		  g_sim = 0;
	 	  tuple_l.clear();
		  //collect all "heavy" neighbor pairs in tuples 
		  forall_adj_nodes(neigh1, v1)
		      forall_adj_nodes(neigh2, v2)
			  if (prev_matrix(index(neigh1), index(neigh2)) > cut_off){
			      n_t.first() = index(neigh1);
			      n_t.second() = index(neigh2);
			      n_t.third() = prev_matrix(index(neigh1), index(neigh2));
			      tuple_l.append(n_t);
			  }
		  if (greedy_choice == TRUE){
		      //greedy set packing for computing local match score
		      tuple_l.sort(comp_triple);
    		      l_m1.clear(); l_m2.clear(); 	  		
		      forall(n_t, tuple_l){
			  l_m1.insert(n_t.first(), FALSE);
			  l_m2.insert(n_t.second(), FALSE);
		      }
		      forall(n_t, tuple_l)
			  if ((l_m1.access(n_t.first()) == FALSE) &&
			      (l_m2.access(n_t.second()) == FALSE)){
			       l_m1.insert(n_t.first(), TRUE);
			       l_m2.insert(n_t.second(), TRUE);
		               g_sim += n_t.third() / 
			                (ppi1.degree(in_to_node1[n_t.first()]) * 
			                 ppi2.degree(in_to_node2[n_t.second()]));
			  }
		      g_sim *= sqrt(tuple_l.length());
	          }//end of local greedy set packing
		  else{
		      local_bg.clear();
		      part1.clear();
		      part2.clear();
		      match_list.clear();
		      from_ppi_to_part1.clear();
		      from_part_to_ppi1.clear();
		      from_ppi_to_part2.clear();
		      from_part_to_ppi2.clear();
		      weight.init(local_bg, tuple_l.length(), 0);
		      //all "heavy" neighbor pairs in local_bg
		      forall(n_t, tuple_l){
			  if (from_ppi_to_part1.defined(in_to_node1[n_t.first()]) == FALSE){
			      temp1 = local_bg.new_node();
		 	      from_ppi_to_part1[in_to_node1[n_t.first()]] = temp1;
			      from_part_to_ppi1[temp1] = in_to_node1[n_t.first()];
			      part1.append(temp1);
	    		  }
	    		  else
			      temp1 = from_ppi_to_part1[in_to_node1[n_t.first()]];
		  	  if (from_ppi_to_part2.defined(in_to_node2[n_t.second()]) == FALSE){
			      temp2 = local_bg.new_node();
			      from_ppi_to_part2[in_to_node2[n_t.second()]] = temp2;
			      from_part_to_ppi2[temp2] = in_to_node2[n_t.second()];
			      part2.append(temp2);
			  }
	    		  else
			      temp2 = from_ppi_to_part2[in_to_node2[n_t.second()]];
			  e = local_bg.new_edge(temp1, temp2);
			  weight[e] = n_t.third();
		      }
    		      MWBM_SCALE_WEIGHTS(local_bg, weight);
		      match_list = MAX_WEIGHT_BIPARTITE_MATCHING(local_bg,part1,part2,weight);
		      forall(e, match_list)
			  g_sim += weight[e] / 
			           (ppi1.degree(from_part_to_ppi1[local_bg.source(e)]) *
			            ppi2.degree(from_part_to_ppi2[local_bg.target(e)]));
		      g_sim *= sqrt(tuple_l.length());	
		  }//end of local optimum set packing
		  (*cur_matrix_p)(in1, in2) = alpha_cons*g_sim + (1- alpha_cons) * (*simmatrix_p)(in1, in2);
	     }//end of all ppi2 nodes
	}//end of all ppi1 nodes, end of local set pack

	//2- update prev_matrix with cur_matrix values
        for (in1=0; in1<len1; in1++)
	    for (in2=0; in2<len2; in2++)
		prev_matrix(in1,in2) = (*cur_matrix_p)(in1, in2);

    }//end of iterations
}//end of optimize_local_score

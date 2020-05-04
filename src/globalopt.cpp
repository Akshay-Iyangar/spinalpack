#include "globalopt.h"

//globals
array2<double>*  cur_mat_p_global;
int len1_global, len2_global;
//hashing related globals and functions
typedef int (*hash_ftype)(int, int, int, int);
typedef two_tuple<int,int> (*rev_hash_ftype)(int, int, int);
hash_ftype hash_f;
rev_hash_ftype rev_hash_f;

int hash1(int in1, int in2, int len1, int len2){
	return (in1*len2 + in2);
}
int hash2(int in1, int in2, int len1, int len2){
	return (in2*len1 + in1);
}
two_tuple<int, int> rev_hash1(int val, int len1, int len2){
	two_tuple<int, int> n_t;
	n_t.first() = val / len2;
	n_t.second() = val - n_t.first()*len2;
	return (n_t);
}
two_tuple<int, int> rev_hash2(int val, int len1, int len2){
	two_tuple<int, int> n_t;
	n_t.second() = val / len1;
	n_t.first() = val - n_t.second()*len1;
	return (n_t);
}

//used for (reverse) sorting. useful in global pack. 
int comp_node_triple_global(const int& hash_val1, const int& hash_val2){
    two_tuple<int, int> n_t1, n_t2;	
    n_t1 = rev_hash_f(hash_val1, len1_global, len2_global);
    n_t2 = rev_hash_f(hash_val2, len1_global, len2_global);

    if ( (*cur_mat_p_global)(n_t1.first(), n_t1.second()) > (*cur_mat_p_global)(n_t2.first(), n_t2.second()) ) 
	return -1;
    else if ( (*cur_mat_p_global)(n_t1.first(), n_t1.second()) < (*cur_mat_p_global)(n_t2.first(), n_t2.second()) ) 
	return 1;
    else 
	return 0;
}

void evaluate_score(h_array<node, node> align_map, graph ppi1, graph ppi2, array2<double>* simmatrix_p, 
		    double alpha_cons, FILE *outfile){
    int shared_ints = 0;
    double seq_sim = 0, score;
    node v1, adj_v1; 
    char match_line[100];

    //evaluate score
    forall_defined(v1, align_map)
	forall(adj_v1, ppi1.adj_nodes(v1))
	    if ((align_map.defined(adj_v1)) && (ppi2.adj_nodes(align_map[v1]).search(align_map[adj_v1]) != nil))
		shared_ints ++;
    forall_defined(v1, align_map)
	seq_sim += (*simmatrix_p)(index(v1), index(align_map[v1]));
    printf("Shared Interactions vs Total Sequence Similarity: %d %f\n", shared_ints/2, seq_sim);
    score = alpha_cons * (shared_ints/2) + (1-alpha_cons) * seq_sim;
    printf("Global Score: %.12f\n", score);

    sprintf(match_line, "!Alpha: %f", alpha_cons);
    fputs(match_line, outfile);
    sprintf(match_line, "!Shared Interactions vs Total Sequence Similarity: %d %f\n", shared_ints/2, seq_sim);
    fputs(match_line, outfile);
    sprintf(match_line, "!Global Score: %.12f\n", score);
    fputs(match_line, outfile);
}

void optimize_global_score2(array2<double>* cur_mat_p, array2<double>* simmatrix_p, graph ppi1, graph ppi2, 
			    double alpha_cons, char *matchfile_name){

    graph bg;
    list<node> part1, part2;
    map<node, node> map1, map1rev, map2, map2rev;
    edge_array<double> weight;
    h_array<node, node> align_map;
    //insert bg nodes wrt ppi1,2 nodes
    node v1, v2, temp;
    edge e;
    forall_nodes(v1, ppi1){
        temp = bg.new_node();
	part1.append(temp);
	map1[v1] = temp;
	map1rev[temp] = v1;
    }
    forall_nodes(v2, ppi2){
        temp = bg.new_node();
	part2.append(temp);
	map2[v2] = temp;
	map2rev[temp] = v2;
    }
    //insert edges in bg and attach weights
    weight.init(bg, ppi1.number_of_nodes()*ppi2.number_of_nodes(), 0);  
    forall_nodes(v1, ppi1)
	forall_nodes(v2, ppi2)
	    if ((*cur_mat_p)(index(v1), index(v2)) > 0)
		{
	    	    e = bg.new_edge(map1[v1], map2[v2]);
	    	    weight[e] = (*cur_mat_p)(index(v1), index(v2));
		}	
    //maxweight bipartite matching
    MWBM_SCALE_WEIGHTS(bg,weight);
    list<edge> M = MAX_WEIGHT_BIPARTITE_MATCHING(bg,part1,part2,weight);
    //insert updated M edges into alignment
    forall(e, M)
	align_map[map1rev[bg.source(e)]] = map2rev[bg.target(e)];
 
    //write alignment output into file
    FILE *outfile;
    outfile = fopen(matchfile_name, "w");
    evaluate_score(align_map, ppi1, ppi2, simmatrix_p, alpha_cons, outfile);

    char match_line[30];
    forall_defined(v1, align_map){
	sprintf(match_line, "%d %d\n", index(v1), index(align_map[v1]));
	fputs(match_line, outfile);
    }
    fclose(outfile);

}
void optimize_global_score(array2<double>* cur_mat_p, array2<double>* simmatrix_p, graph ppi1, graph ppi2, 
			    double alpha_cons, char *matchfile_name){
    //assign globals
    len1_global = ppi1.number_of_nodes();
    len2_global = ppi2.number_of_nodes();
    if (len1_global > len2_global) {hash_f = &hash1; rev_hash_f = &rev_hash1;}
    else {hash_f=&hash2; rev_hash_f = &rev_hash2;}
    cur_mat_p_global = cur_mat_p;

    edge e, bg_e, alt_e1, alt_e2;
    graph bg;
    node v1, v2, best1, best2, adj_v1, adj_v2, temp;
    map<node, bool> aligned1, aligned2;
    h_array<node, node> align_map;
    list<edge> M;
    map<edge, bool> edge_in_M;
    dictionary<node, edge> node_to_M;
    three_tuple<node, node, double> edge_tup;
    two_tuple<int, int> index_tup;
    list<int> unaligned_pairs;
    map<int, node> in_to_node1, in_to_node2;

    list<node> part1, part2, prev1, prev2;
    map<node, node> bg_to_ppi1, bg_to_ppi2, ppi1_to_bg, ppi2_to_bg;
    edge_array<double> weight;
    map2<node, node, int> freq;
    double max_val, e_gain, overlap_gain; 

    //index-to-node map
    forall_nodes(v1, ppi1)
	in_to_node1[index(v1)] = v1;
    forall_nodes(v2, ppi2)
	in_to_node2[index(v2)] = v2;
    //construct unaligned_pairs list and sort wrt similarity values computed in local_optimization
    forall_nodes(v1, ppi1)
	forall_nodes(v2, ppi2){
	    unaligned_pairs.append(hash_f(index(v1), index(v2), len1_global, len2_global));
	}
    unaligned_pairs.sort(comp_node_triple_global);
    //iteratively find alignment seed pair and expand alignment
    while (unaligned_pairs.empty() == FALSE){
	//find initial pair with the heaviest similarity among all unmatched_pairs  
	while(unaligned_pairs.empty() == FALSE){
	    index_tup = rev_hash_f(unaligned_pairs.pop(), len1_global, len2_global);
	    edge_tup.first() = in_to_node1[index_tup.first()];
	    edge_tup.second() = in_to_node2[index_tup.second()];
 	    edge_tup.third() = (*cur_mat_p)(index_tup.first(), index_tup.second());
	    if ( (aligned1[edge_tup.first()] == FALSE) && (aligned2[edge_tup.second()] == FALSE) ){
		best1 = edge_tup.first();
		best2 = edge_tup.second();
		max_val = edge_tup.third();
		break;
	    }
	}
	if ((unaligned_pairs.empty()) || (max_val == 0))
	    break;
	//insert heaviest pair into alignment
        align_map[best1] = best2;
        aligned1[best1] = TRUE;
        aligned2[best2] = TRUE;
        temp = bg.new_node();
        part1.append(temp); bg_to_ppi1[temp] = best1;
        temp = bg.new_node();
        part2.append(temp); bg_to_ppi2[temp] = best2;
	//expand alignment each time considering neighbors of aligned nodes until none
        while ((part1.empty() == FALSE) && (part2.empty() == FALSE)){
	    prev1.clear();
  	    forall(v1, part1)
	        prev1.append(bg_to_ppi1[v1]);
	    prev2.clear();
	    forall(v2, part2)
	        prev2.append(bg_to_ppi2[v2]);

	    bg.clear();
	    part1.clear();
	    part2.clear();
	    bg_to_ppi1.clear();
	    bg_to_ppi2.clear();
	    ppi1_to_bg.clear();
	    ppi2_to_bg.clear();
	    freq.clear();
	    //create node sets in bg: part1, part2 and mappings bw bg nodes and ppi nodes
	    forall(v1, prev1){
		//if v1 from previous bg was not aligned add it to this bg as well
	        if ((aligned1[v1] == FALSE) && (ppi1_to_bg.defined(v1) == FALSE)){
        	    temp = bg.new_node();
		    part1.append(temp);
		    bg_to_ppi1[temp] = v1;
		    ppi1_to_bg[v1] = temp;
	        }
		//add all neighbors of v1 into bg if they are not already aligned
	        forall(adj_v1, ppi1.adj_nodes(v1))
		    if ((aligned1[adj_v1] == FALSE) && (ppi1_to_bg.defined(adj_v1) == FALSE)){
        	        temp = bg.new_node();
		        part1.append(temp);
		        bg_to_ppi1[temp] = adj_v1;
		        ppi1_to_bg[adj_v1] = temp;
	            }
	    }//end of v1 processing
	    forall(v2, prev2){
		//if v2 from previous bg was not aligned add it to this bg as well
	        if ((aligned2[v2] == FALSE) && (ppi2_to_bg.defined(v2) == FALSE)){
        	    temp = bg.new_node();
		    part2.append(temp);
		    bg_to_ppi2[temp] = v2;
		    ppi2_to_bg[v2] = temp;
	        }
		//add all neighbors of v2 into bg if they are not already aligned
	        forall(adj_v2, ppi2.adj_nodes(v2))
		    if ((aligned2[adj_v2] == FALSE) && (ppi2_to_bg.defined(adj_v2) == FALSE)){
        	        temp = bg.new_node();
		        part2.append(temp);
		        bg_to_ppi2[temp] = adj_v2;
		        ppi2_to_bg[adj_v2] = temp;
	            } 
	    }//end of v2 processing

	    weight.init(bg, part1.length()*part2.length(), 0); 
	    //insert edges in bg
	    forall(v1, part1)
	        forall(adj_v1, ppi1.adj_nodes(bg_to_ppi1[v1]))
		    if (align_map.defined(adj_v1))
		        forall(v2, ppi2.adj_nodes(align_map[adj_v1]))
			    if (part2.search(ppi2_to_bg[v2]) != nil)
				//if occured before, increment frequency
				//freq denotes the number of shared interactions
				//of v1, v2 with the previously aligned nodes
				if (freq.defined(v1, ppi2_to_bg[v2]))
				    freq(v1, ppi2_to_bg[v2]) ++;
				//if first time, insert into bg
				else{
			            e = bg.new_edge(v1, ppi2_to_bg[v2]);
			            weight[e] = (*cur_mat_p)(index(bg_to_ppi1[v1]), index(v2));
				    freq(v1, ppi2_to_bg[v2]) = 1;
				}

    	    //maxweight bipartite matching
	    MWBM_SCALE_WEIGHTS(bg,weight);
	    M.clear();
            M = MAX_WEIGHT_BIPARTITE_MATCHING(bg, part1, part2, weight);
	    if (M.length() == 0)
		break;
	    //improve matching to be added into alignment
	    node_to_M.clear();
	    forall_edges(e, bg)
		edge_in_M[e] = FALSE;
	    forall(e, M){
		node_to_M.insert(bg.source(e), e);
		node_to_M.insert(bg.target(e), e);
		edge_in_M[e] = TRUE;
	    }
	    forall_edges(e, bg)
		if (edge_in_M[e] == FALSE)
		    //if edge e overlaps with two edges in M
		    if ( (node_to_M.lookup(bg.source(e)) != nil) && (node_to_M.lookup(bg.target(e)) != nil) ){
			alt_e1 = node_to_M[node_to_M.lookup(bg.source(e))]; 
			alt_e2 = node_to_M[node_to_M.lookup(bg.target(e))];
			//compute gain of replacing overlaps of e with e
			if (alpha_cons >= 0.5){
			    e_gain = freq(bg.source(e), bg.target(e));
			    overlap_gain = freq(bg.source(alt_e1), bg.target(alt_e1)) + freq(bg.source(alt_e2), bg.target(alt_e2));
			}
			else{
			    e_gain = alpha_cons * freq(bg.source(e), bg.target(e)) + 
				    (1-alpha_cons) * (*simmatrix_p)(index(bg_to_ppi1[bg.source(e)]), index(bg_to_ppi2[bg.target(e)]));
			    overlap_gain = alpha_cons * 
					   (freq(bg.source(alt_e1), bg.target(alt_e1)) + freq(bg.source(alt_e2), bg.target(alt_e2))) + 
					   (1-alpha_cons) * 
					   ((*simmatrix_p)(index(bg_to_ppi1[bg.source(alt_e1)]), index(bg_to_ppi2[bg.target(alt_e1)])) + 
					    (*simmatrix_p)(index(bg_to_ppi1[bg.source(alt_e2)]), index(bg_to_ppi2[bg.target(alt_e2)])));
			}
			//replacing overlaps of e in M, with e increases alignment score
		        if ( e_gain >= overlap_gain ){
	        	    node_to_M.del_item(node_to_M.lookup(bg.source(alt_e1)));
	        	    node_to_M.del_item(node_to_M.lookup(bg.target(alt_e1)));
			    edge_in_M[alt_e1] = FALSE;
			    M.remove(alt_e1);
	        	    node_to_M.del_item(node_to_M.lookup(bg.source(alt_e2)));
	        	    node_to_M.del_item(node_to_M.lookup(bg.target(alt_e2)));
			    edge_in_M[alt_e2] = FALSE;
			    M.remove(alt_e2);
			    node_to_M.insert(bg.source(e), e);
			    node_to_M.insert(bg.target(e), e);
			    edge_in_M[e] = TRUE;
			    M.append(e);
			}
		    }
		    //if edge e overlaps with only one edge in M
		    else if ((node_to_M.lookup(bg.source(e)) != nil) || (node_to_M.lookup(bg.target(e)) != nil)){
			if (node_to_M.lookup(bg.source(e)) != nil) 
			    alt_e1 = node_to_M[node_to_M.lookup(bg.source(e))]; 
			else
			    alt_e1 = node_to_M[node_to_M.lookup(bg.target(e))];
			//compute gain of replacing single overlap of e with e
			if (alpha_cons >= 0.5){
			    e_gain = freq(bg.source(e), bg.target(e));
			    overlap_gain = freq(bg.source(alt_e1), bg.target(alt_e1));
			}
			else{
			    e_gain = alpha_cons * freq(bg.source(e), bg.target(e)) + 
				    (1-alpha_cons) * (*simmatrix_p)(index(bg_to_ppi1[bg.source(e)]), index(bg_to_ppi2[bg.target(e)]));
			    overlap_gain = alpha_cons * freq(bg.source(alt_e1), bg.target(alt_e1)) + 
					   (1-alpha_cons) * 
					   (*simmatrix_p)(index(bg_to_ppi1[bg.source(alt_e1)]), index(bg_to_ppi2[bg.target(alt_e1)]));
			}
			//replacing single overlap of e in M, with e increases alignment score
		        if ( e_gain > overlap_gain ){
	        	    node_to_M.del_item(node_to_M.lookup(bg.source(alt_e1)));
	        	    node_to_M.del_item(node_to_M.lookup(bg.target(alt_e1)));
			    edge_in_M[alt_e1] = FALSE;
			    M.remove(alt_e1);
			    node_to_M.insert(bg.source(e), e);
			    node_to_M.insert(bg.target(e), e);
			    edge_in_M[e] = TRUE;
			    M.append(e);
			}
		    }
		    //if edge e overlaps with no edge in M
		    //do not automatically insert in M, it will be 
		    //considered with the next batch

	    //insert updated M edges into alignment
	    forall(e, M){
	    	best1 = bg_to_ppi1[bg.source(e)];
	        best2 = bg_to_ppi2[bg.target(e)];
	        align_map[best1] = best2;
    	        aligned1[best1] = TRUE;
    	        aligned2[best2] = TRUE;
	    }
        }//end while first seed

	//clear everything and continue next iteration
	bg.clear();
	part1.clear();
	part2.clear();
	bg_to_ppi1.clear();
	bg_to_ppi2.clear();
	ppi1_to_bg.clear();
	ppi2_to_bg.clear();
	freq.clear();
    }//end while all seeds

    //write alignment output into file
    FILE *outfile;
    outfile = fopen(matchfile_name, "w");
    evaluate_score(align_map, ppi1, ppi2, simmatrix_p, alpha_cons, outfile);
    char match_line[30];
    forall_defined(v1, align_map){
	sprintf(match_line, "%d %d\n", index(v1), index(align_map[v1]));
	fputs(match_line, outfile);
    }
    fclose(outfile);

}


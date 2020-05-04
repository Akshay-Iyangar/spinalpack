#include "common.h"
#include "initialize.h"
#include "localopt.h"
#include "globalopt.h"

int main(int argc, char* argv[])
{
    //compile with:
    //g++-4.1 -I$LEDAROOT/incl -L$LEDAROOT initialize.cpp localopt.cpp globalopt.cpp ppialign.cpp -lGeoW -lD3 -lW -lP -lG -lL -lX11 -lm -O1 -o spinal
    //sample run with seq_sim matrix and default parameters: 
    //./spinal -II -ns ./data/dmela.tab.gml ./data/scere.tab.gml ./data/dmela-scere.evals.pin ./data/results.txt 0.7 
    //sample run with no seq_sim matrix and default parameters:
    //./spinal -II -n ./data/dmela.tab.gml ./data/scere.tab.gml ./data/results.txt  
    char *ppi1file_name, *ppi2file_name, *seqsimfile_name, *matchfile_name, *align_choice, *global_choice;
    double alpha_cons;
    bool local_pack_ver;
    int local_iter_cnt, local_cut_offpnt;

    global_choice = argv[1];
    align_choice = argv[2];
    ppi1file_name = argv[3];
    ppi2file_name = argv[4];

    //create ppis, remove selfloops and make undirected	
    printf("Reading PPI Networks...\n");
    graph ppi1, ppi2;
    ppi1.read_gml(ppi1file_name);
    ppi2.read_gml(ppi2file_name);
    ppi1.make_undirected();
    ppi2.make_undirected();
    int len1 = ppi1.number_of_nodes();
    int len2 = ppi2.number_of_nodes();

    if (strcmp(align_choice, "-n") == 0){
	matchfile_name = argv[5];
	alpha_cons = 1;
	if (argc > 6){
            //TRUE greedy matching, FALSE optimum matching
            local_pack_ver = atoi(argv[6]);
            local_iter_cnt = atoi(argv[7]);
            local_cut_offpnt = atoi(argv[8]);
	}
	else{
            local_pack_ver = TRUE;
            local_iter_cnt = 14;
	    if (len1*len2 <= 200000)
		local_cut_offpnt = (len1*len2)/2;
	    else
                local_cut_offpnt = 100000;
	}
    }
    else if (strcmp(align_choice, "-ns") == 0){
	seqsimfile_name = argv[5];
	matchfile_name = argv[6];
	alpha_cons = atof(argv[7]);
	if (argc > 8){
            //TRUE greedy matching, FALSE optimum matching
            local_pack_ver = atoi(argv[8]);
            local_iter_cnt = atoi(argv[9]);
            local_cut_offpnt = atoi(argv[10]);
	}
	else{
            local_pack_ver = TRUE;
            local_iter_cnt = 14;
	    if (len1*len2 <= 200000)
		local_cut_offpnt = (len1*len2)/2;
	    else
                local_cut_offpnt = 100000;
	}
    }
    else{
	printf("Alignment choice must be one of -n or -ns");
	exit(1);
    }
 

    //construct and process sequence similarity matrix
    array2<double> simmatrix(len1, len2);  
    if (strcmp(align_choice, "-ns") == 0){
        printf("Constructing Sequence Similarity Matrix...\n");
	construct_simmatrix(&simmatrix, seqsimfile_name);
    }

    //initialize current matrix of similarity(wrt objective func.) values
    printf("Initializing Similarity Matrix...\n");
    array2<double> cur(len1, len2);
    init_cur(&cur, &simmatrix, ppi1, ppi2, alpha_cons);


    //optimize local score function
    printf("Starting Alignment:\n");
    optimize_local_score(&cur, &simmatrix, ppi1, ppi2, alpha_cons, 
			 local_iter_cnt, local_cut_offpnt, local_pack_ver);

    //optimize global score function
    printf("Wrapping Up...\n");
    if (strcmp(global_choice, "-I") == 0)
        optimize_global_score2(&cur, &simmatrix, ppi1, ppi2, alpha_cons, matchfile_name); 
    else if (strcmp(global_choice, "-II") == 0)
        optimize_global_score(&cur, &simmatrix, ppi1, ppi2, alpha_cons, matchfile_name); 

}

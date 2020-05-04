#include "initialize.h"

//construct similarity matrix from file
void construct_simmatrix(array2<double>* sim_matrix_p, char* seqsimfile_name){
    char simline[60];
    int n1, n2;
    double simval; 
    FILE *simfile; 
    simfile = fopen(seqsimfile_name, "r");

    while(fgets(simline,sizeof(simline),simfile) != NULL){
	// strip trailing '\n' 
	if(simline[strlen(simline)-1] == '\n') 
	    simline[strlen(simline)-1] = 0;
	n1 = atoi(strtok(simline, " "));
	n2 = atoi(strtok(NULL, " "));
	simval = strtod(strtok(NULL, " "), NULL);
	//assign similarity matrix value
	(*sim_matrix_p)(n1, n2) = simval;
    }
    fclose(simfile);
}

//initialize cur matrix
void init_cur(array2<double>* cur_matrix_p, array2<double>* sim_matrix_p, graph ppi1, graph ppi2, double alpha_cons){ 
    node v1, v2;
    int maxdeg1 = 0, degsum1 = 0, avgdeg1 = 0;
    int mindeg1 = ppi1.degree(ppi1.first_node());
    forall_nodes(v1, ppi1){
	degsum1 += ppi1.degree(v1);
	if (ppi1.degree(v1) > maxdeg1) maxdeg1 = ppi1.degree(v1);
	if (ppi1.degree(v1) < mindeg1) mindeg1 = ppi1.degree(v1);
    }
    avgdeg1 = degsum1 / ppi1.number_of_nodes(); 

    int maxdeg2 = 0, degsum2 = 0, avgdeg2 = 0;
    int mindeg2 = ppi2.degree(ppi2.first_node());
    forall_nodes(v2, ppi2){
	degsum2 += ppi2.degree(v2);
	if (ppi2.degree(v2) > maxdeg2) maxdeg2 = ppi2.degree(v2);
	if (ppi2.degree(v2) < mindeg2) mindeg2 = ppi2.degree(v2);
    }
    avgdeg2 = degsum2 / ppi2.number_of_nodes();

    double avgdegratio;
    if (avgdeg1 > avgdeg2) avgdegratio = avgdeg1 / avgdeg2;
    else avgdegratio = avgdeg2 / avgdeg1;

    int max_deg_diff = max(maxdeg1 - mindeg2, maxdeg2 - mindeg1);
    double deg_diff; 
    forall_nodes(v1, ppi1)
	forall_nodes(v2, ppi2){
	    deg_diff = 1-(abs(ppi1.degree(v1) - ppi2.degree(v2)) / max_deg_diff);
	    if (alpha_cons == 1)
		deg_diff *= (ppi1.degree(v1) + ppi2.degree(v2)) / (maxdeg1 + maxdeg2); 	
	    (*cur_matrix_p)(index(v1), index(v2)) = alpha_cons*deg_diff + 
			                            (1- alpha_cons) * (*sim_matrix_p)(index(v1), index(v2));
    	}
}

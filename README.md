SPINAL-Scalable PPI Network Alignment, by C. Erten.

SPINAL has been tested on Unix systems. It has been coded in C++ and compiled via GNU compiler, g++4.1. We provide several useful Python scripts as well. You may need the “NetworkX” package for some of the Python scripts. Details of this package can be found at: http://networkx.lanl.gov/

We note that checking the format of the files in the data folder after reading the execution instructions below might be quite helpful.

USAGE:
- If using a sequence similarity file, SPINAL can be run with default parameters 
  as follows:

  ./spinal -ver -ns X.gml Y.gml seq_file out_file alpha

  where -ver is either “-I” (if only running Algorithm1) or “-II” (if running   
  both Algorithm1 and 2), X.gml is the ppi network of organism1, Y.gml is that 
  of organism2, seq_file is the normalized sequence similarity scores, out_file 
  is the output file where the alignment results are written. 

  NOTE1: X.gml, Y.gml must be in .gml format. In this format, each node of X.gml 
  must have integer indices starting with 0, and “labels” corresponding to  
  protein names. File seq_file is assumed to be consisting of lines 
  “x1 y1   score” where x1 is the “index” of a protein in X.gml and x2 is the 
  “index” of a protein in Y.gml, score is a normalized sequence similarity score 
  between 0 and 1. To help for construction of X.gml, Y.gml, seq_file in the 
  desired format we provide a python script, data_generation.py. If you already 
  have a ppi file as a text file in “exactly” the following format:
	INTERACTOR_A  INTERACTOR_B
	x1 x2
	x1 x5
	 ... 
  where x1, x2, etc. are protein ids of your choice, and a sequence similarity 
  file in the following format:
	x1 y1 40
	x1 y5 55
   	  ...
  where x1, y1 etc. are protein ids as before and the numbers are similarity 
  values, you can run, 
  ./pyscripts/data_generation.py X.txt Y.txt S
  where X.txt, Y.txt are ppi networks in the defined format and S is the 
  sequence similarity file. This generates X.gml, Y.gml, and S.pin which are  
  files in the desired format. You can then use these as input to SPINAL.

  NOTE2: The output of SPINAL, “out_file” provides alignments, but with indices 
  in the corresponding .gml files. If you want to convert back to original names 
  in X.txt and Y.txt, you may use the spinal_original_names.py as follows:
  .pyscripts/spinal_original_names.py X.gml Y.gml out_file out_file_original  

- If you do not want to use a sequence similarity file, default parameter 
  execution is as follows:
  ./spinal -ver -n p1 p2 out_file alpha
  The discussion above applies similarly. 

SAMPLE RUNS:
- Sample run with sequence similarity file and default parameters: 
  ./spinal -II -ns ./data/dmela/dmela.tab.gml ./data/scere/scere.tab.gml  \     
           ./data/dm-sc/dmelascere.evals.pin ./data/dm-sc/results.txt 0.7 

- Sample run with no sequence similarity file and default parameters:
  ./spinal -II -n ./data/dmela/dmela.tab.gml ./data/scere/scere.tab.gml   \
  			./data/dm-sc/results.txt 


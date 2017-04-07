# Simple_Gfa_parser

Script used in order to parse assembly graph in gfa format and store it in Networkx and graph-tool datastructure.
Depending on entry argument, it can also generate a gfa file from the fasta file of contigs from assembly, in order to do this it call in turn megahit and Bandage which need to by availlable in your Path. Graph object are stored as the native .gt format for graphtool and as a pickle output for networkx. 

An important point of this script is that we delete a few edges in order to separate contigs from their reverse complements. This is done by looking at the shortest paht between a contig and its reverse complement and deleting the shortest among all contigs. 

Exemple of use : Parse_gfa.py overlaping_size file_in file_out -C contig_to_bin_map
    overlaping_size :  number of nucleotides overlaping between contigs in the assembly graph
    file_in : if file_in is ".fa", call megahit toolbox then Bandage to produce a gfa file, then process it
              if file_in is ".fastg", call Bandage to to produce a gfa file, then process it
              if file_in is ".gfa", just extract graph informations
    file_out : name of output file
    -C contig_to_bin_map : allow to add information regarding vertexes in your graph, namely membership of a contig in a bin. Need to be in the following format : name_of_contig:name_of_bin .Name_of_contig need to be the same as in file_in. 
    


# Simple_Gfa_parser

This script allow to parse assembly graph in .fa, .fastg and .gfa format, it is then stored in Networkx and graph-tool datastructure.
Depending on entry argument, gfa file is generated from .fa or .fastg. In order to do so,  megahit and Bandage are called in turn. This script thus requiere them to by availlable in your Path. Graph objects are stored as the native .gt format for graphtool and as a pickle output for networkx. 

An important feature of this script is that, if needed, a few edges are deleted in order to separate contigs from their reverse complements. This is done by looking at the collection of shortest path between a contig and their reverse complement, then the shortest path of the collection is cut at the vertex of lowest degre. 

Exemple of use : 
```
Parse_gfa.py overlaping_size file_in file_out -C contig_to_bin_map -Draw Graph.format
```
overlaping_size :  number of nucleotides overlaping between contigs in the assembly graph

file_in : 
* if file_in is ".fa", call megahit toolbox then Bandage to produce a gfa file, then process it
* if file_in is ".fastg", call Bandage to to produce a gfa file, then process it
* if file_in is ".gfa", just extract graph informations

file_out : name of output file

-C contig_to_bin_map : falcultative argument, allow to add information regarding vertexes in your graph, namely membership of a contig in a bin. The file contig_to_bin_map need to be in the following format : name_of_contig:name_of_bin .Name_of_contig need to be the same as in file_in. At the moment a maximum of 20 differents bin is supported. 
    
-Draw Graph.format : facultative argument, use graph-tool to draw the assembly graph in the file "graph" in the format ".format", it can be any of thoses : xlib", "ps", "svg", "svgz", "fig", "mif", "hpgl", "pcl", "png", "gif", "dia", "imap", "cmapx". If contig_to_bin_map was also used, vertices are colored following vertices bin membership. 

#Pratical example :
Using files in Example folder, the following command is ewecuted :
```
Parse_gfa.py 141 Example/LactoParaPhage_contig.fa Example/LactoParaPhage_contig -C Lacto_Species.csv -Draw Graph.pdf
```
This yield the following graph-tool graph :
![picture alt](https://github.com/Sebastien-Raguideau/Simple_Gfa_parser/blob/master/Example/Graph.pdf)

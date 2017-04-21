#!/usr/bin/env python
# -*- coding: latin-1 -*-
import os 
import networkx as nx
import Bio
from Bio.SeqIO.FastaIO import *
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import argparse
import numpy as np
from collections import Counter
import graph_tool as gt
from graph_tool import draw
import pickle

def Parse_gfa(gfa_file,G,Dictionary_contig_label):
	def reverse_sign(sign):
		sign="+"*(sign=="-")+"-"*(sign=="+")
		return sign 
	Handle=open(gfa_file)
	for Line in Handle : 
		Line=Line.rstrip().split("\t")
		if Line[0]=="S" : 
			# [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
			(name,seq_plus)=(Line[1],Seq(Line[2]))
			seq_len=len(seq_plus)
			color="#CF0626"
			if len(Line)>5 :
				color=Line[5].split(':')[-1]
			if name in Dictionary_contig_label :
				Bin=Dictionary_contig_label[name]
			else :
				Bin=[]
			cov=float(Line[4].replace("KC:i:",""))/seq_len
			G.add_node(name+"+",Name=name+"+",Seq=seq_plus,Cov=cov,Color=color,Bin=Bin)
			G.add_node(name+"-",Name=name+"-",Seq=seq_plus.reverse_complement(),Cov=cov,Color=color,Bin=Bin)
		if Line[0]=="L" :
			# ['L','2','+','27','+','99M']
			G.add_edge("".join(Line[1:3]),"".join(Line[3:5]))
			G["".join(Line[1:3])]["".join(Line[3:5])]['Struct']="\t".join(Line)
			G.add_edge("".join([Line[3],reverse_sign(Line[4])]),"".join([Line[1],reverse_sign(Line[2])]))
	Handle.close()

def Generate_gfa_file(G,gfa_file,contig_overlap,After_pruning) :
	Handle=open(gfa_file,"w")
	for node_name,dico_node in G.node.items() :
	#S	1	AGTCTTCGTCCAGGGGGCCGCCTTCGCCACCGGTATTCCTCCAGATCTCTACGCATTTCACCGCTACACCTGGAATTCTACCCCCCTCTACGAGACTCACGCTTGCCAGTATCAGATG	LN:i:118	KC:i:82948	CL:z:#5cb9be	C2:z:#5cb9be
		if node_name[-1]=="+" :
			seq=dico_node["Seq"]
			color=""
			if 'Color' in dico_node :
				color="\tCL:z:"+dico_node['Color']
			towrite=["S",node_name[:-1],seq.tostring(),"LN:i:"+str(len(seq)),"KC:i:"+str(len(seq)*dico_node['Cov'])]
			Handle.write("\t".join(towrite)+color+"\n")
	for Random_Edges in G.edges() :
	#L	52	-	11	+	99M
		if After_pruning:
			node1,node2=(G.node[Random_Edges[0]],G.node[Random_Edges[1]])
			if node1['Seq'][-int(contig_overlap):]==node2['Seq'][:int(contig_overlap)] :
				Edges=[node1['Name'],node2['Name']]
			else :
				Edges=[node2['Name'],node1['Name']]
			towrite=["L",Edges[0][:-1],Edges[0][-1],Edges[1][:-1],Edges[1][-1],contig_overlap+"M\n"]
			Handle.write("\t".join(towrite))
		else :
			Dico=G[Random_Edges[0]][Random_Edges[1]]
			if Dico!={} :
				Handle.write(Dico["Struct"]+"\n")
	Handle.close()


# def Generate_adjacency_matrix(G) :
# 	Mat=np.zeros((len(G),len(G)))
# 	Ordered_node_list=sorted(G.node.keys(),key=lambda x:int(x.split("_")[-1][:-1]))
# 	for index,node in enumerate(Ordered_node_list) :
# 		for neighbors in G.neighbors(node) :
# 			Mat[index][Ordered_node_list.index(neighbors)]=1
# 	return Mat,Ordered_node_list

# def Print_adjacency_matrix(G,file) :
# 	Mat,Ordered_node_list=Generate_adjacency_matrix(G)
# 	Handle=open(file,'w')
# 	Handle.write('\t'.join(Ordered_node_list)+'\n')
# 	Handle.write("\n".join(['\t'.join([str(nb) for nb in line])for line in Mat]))
# 	Handle.close()


def Delete_root_of_joined_comp(G,comp):
	"""For some contigs, paths can be found between themselves and their respective reverse-complement, we make the decision to cut them"""
	def reverse_Node_sign(Node):
		sign=Node[-1]
		Node=Node[:-1]
		sign="+"*(sign=="-")+"-"*(sign=="+")
		Node=Node+sign
		return Node
	def check_result_of_removal(G,comp,Edges_to_remove):
		G_sub=G.subgraph(comp)
		G_sub.remove_edges_from(Edges_to_remove)
		if len([comp for comp in nx.connected_components(G_sub)])!=2:
			return 0,G_sub
		else :
			return 1,G_sub
	def Get_root_of_joined_comp(G,comp) :
		dico_node_selfdist={node[:-1]:0 for node in list(comp)}
		dico_node_selfdist={node:len(nx.shortest_path(G,node+"+",node+"-")) for node in dico_node_selfdist}
		dico_selfdist_node={dist:[] for dist in dico_node_selfdist.values()}
		for node,dist in dico_node_selfdist.items() : 
			dico_selfdist_node[dist].append(node)
		minDist=min(dico_node_selfdist.values())
		List_nodes=dico_selfdist_node[minDist]
		Root=sorted(List_nodes,key=lambda x:len(G.neighbors(x+'+')))[0]+"+"
		return Root
	def Pruning_routine(G,comp,Edges_to_remove,List_Root) :
		#print len(List_Root)
		Current_Root=Get_root_of_joined_comp(G,comp)
		List_Root.append(Current_Root)
		neighbors=[nodes for nodes in G.neighbors(Current_Root)]
		neighbors_num=Counter(List_Root)[Current_Root]-1
		Edges_to_remove.append((Current_Root,neighbors[neighbors_num]))
		Edges_to_remove.append((reverse_Node_sign(Current_Root),reverse_Node_sign(neighbors[neighbors_num])))
		res,G_sub=check_result_of_removal(G,comp,Edges_to_remove)
		if res!=1 :
			Pruning_routine(G_sub,comp,Edges_to_remove,List_Root)
	List_Root=[]
	Edges_to_remove=[]
	Pruning_routine(G,comp,Edges_to_remove,List_Root)
	#return check_result_of_removal(G,comp,Edges_to_remove)
	G.remove_edges_from(Edges_to_remove)


def delete_redundant_comp(G,round) :
	""" The way graph information is stocked imply that we first need to consider every contigs and their reverse complement before being able to maybe delete half of the nodes"""
	def reverse_Node_sign(Node):
		sign=Node[-1]
		Node=Node[:-1]
		sign="+"*(sign=="-")+"-"*(sign=="+")
		Node=Node+sign
		return Node
	def reverse_set_sign(Set):
		Set=set([reverse_Node_sign(Node) for Node in Set])
		return Set
	Node_to_delete=[]
	Comp_to_check=[]
	Listcomp=[comp for comp in nx.connected_components(G)]
	Listlencomp=[len(comp) for comp in Listcomp]
	Dico_len_comp={i:[] for i in Listlencomp}
	for comp in Listcomp :
		Dico_len_comp[len(comp)].append(comp)
	for Len,List_set_comp in Dico_len_comp.items():
		for index1,set_comp in enumerate(List_set_comp) :
			if set_comp!={} :
				if Len==1 :
					comp=list(set_comp)[0][:-1]+"-"
					if comp not in Node_to_delete : 
						Node_to_delete.append(comp)
				else :
					tocheck=[index for index,set_comp2 in enumerate(List_set_comp) if (reverse_set_sign(set_comp)==set_comp2)&(index!=index1)]
					if len(tocheck)==1 :
						set_comp2=List_set_comp[tocheck[0]]
						List_set_comp[tocheck[0]]={}
						for node in list(set_comp2) :
							Node_to_delete.append(node)
					else :
						Comp_to_check.append(set_comp)
	G.remove_nodes_from(Node_to_delete)
	# return Node_to_delete,Comp_to_check
	if round==0 :
		for comp in Comp_to_check :
			Delete_root_of_joined_comp(G,comp)
		delete_redundant_comp(G,1)

def Get_Contig_Assignment(contig_assignment) :
	Dico_Contig_Bug={}
	if contig_assignment[0]=="C" :
		Handle=open(contig_assignment[1])
		for line in Handle : 
			(contig_name,List_bug_name)=line.rstrip().split(',')[0],line.rstrip().split(',')[1:]
			contig_name=contig_name.split('.')[0]
			if contig_name not in Dico_Contig_Bug :
				Dico_Contig_Bug[contig_name]=List_bug_name
			else :
				Dico_Contig_Bug[contig_name].append(List_bug_name)
		Handle.close()
		return Dico_Contig_Bug
	elif contig_assignment[0]=="TabC" :
		Handle=open(contig_assignment[1])
		Header=Handle.next().rstrip().split(',')[1:]
		for line in Handle :
			line=line.rstrip().split(',')
			key=line[0]
			value=[Header[index] for index,value in enumerate(line[1:]) if float(value)!=0]
			value=(value==[])*["NA"]+(value!=[])*value
			Dico_Contig_Bug[key]=value
	return Dico_Contig_Bug



		



def Rewrite_gfa(Dico_Contig_Bug,fasta_file,gfa_file) :
		# get contig name since megahit_toolkit delete them (grr)
	def Get_contig_id(Handle):
		Ordered_Contig_id=[]
		for title, sequence in SimpleFastaParser(Handle) :
			Ordered_Contig_id.append(title.split(' ')[0])
		return Ordered_Contig_id
	Handle=open(fasta_file)
	Dico_Order_Contig={str(Order+1):Name for Order,Name in enumerate(Get_contig_id(Handle))}
	Handle.close()
	# First assign color to contig 
	#	Color_scheme=[Red1,Green1,BLue1,Pink1,Orange1,Yellow1,Gray1,Cyan1,strange,orangish,purple,red2,Green2,blue2,pink2,Gray2,Orange2,Green3,Cyan2,purple2]
	#["#ff0000","#33cc33","#0066ff","#e600ac","#ff9933","#ffff00","#cccccc","#00ffcc","#666699","#802000","#cc0099","#b32400","#00ff00","#003d99","#ff33cc","#737373","#e6b800","#ccff33","#47d1d1","#990099"]

	if Dico_Contig_Bug!={} :
		Color_scheme=["#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00"]
		def merge_color(Listcolor) :
			total_color=np.zeros(3)
			for color in Listcolor :
				total_color=total_color+np.array([int(color[1:3],16),int(color[3:5],16),int(color[5:],16)])
			int_to_hex=lambda x:hex(int(x))[2:].upper() 
			Correct_int_to_hex=lambda x:int_to_hex(x)*(int_to_hex(x)!="0")+"00"*(int_to_hex(x)=="0")
			Merged_color="#"+"".join([Correct_int_to_hex(value) for value in total_color/len(Listcolor)])
			return Merged_color
		List_Bugs=list(set([z for a in Dico_Contig_Bug.values() for z in a ]))
		if "NA" in List_Bugs :
			del List_Bugs[List_Bugs.index("NA")]
		Dico_Contig_color={Contig:merge_color([Color_scheme[List_Bugs.index(Bug)] for Bug in BugList]) if BugList!=["NA"] else "#000000" for Contig,BugList in Dico_Contig_Bug.items()}
		Dico_Order_Color={Order:Dico_Contig_color[Dico_Order_Contig[Order]] for Order,contig_name in Dico_Order_Contig.items()}
		# next add this information to the gfa file 
	Handle=open(gfa_file)
	NewGfa=""
	for line in Handle :
		line=line.rstrip().split('\t')
		if line[0]=="S" :
			if Dico_Contig_Bug!={} : 
				line.append("CL:z:"+Dico_Order_Color[line[1]]+"\n")
			else :
				line[-1]+='\n'
			line[1]=Dico_Order_Contig[line[1]]
		elif line[0]=="L":
			line[1]=Dico_Order_Contig[line[1]]
			line[3]=Dico_Order_Contig[line[3]]
			line.append('\n')
		NewGfa+="\t".join(line)
	Handle.close()
	H=open(gfa_file,'w')
	H.write(NewGfa)
	H.close()

def Translation_from_NX_to_Gt(G,Gt) :
	Ordered_node_list=sorted(G.node.keys(),key=lambda x:int(x.split("_")[-1][:-1]))
	Gt.add_vertex(len(G))
	vertex_prop_Name = Gt.new_vertex_property("string")
	vertex_prop_Seq = Gt.new_vertex_property("object")
	vertex_prop_Cov = Gt.new_vertex_property("double")
	vertex_prop_Color = Gt.new_vertex_property("string")
	vertex_prop_Bin = Gt.new_vertex_property("vector<string>")
	for index,node_name in enumerate(Ordered_node_list) :
		vertex=Gt.vertex(index)
		vertex_prop_Name[vertex]=G.node[node_name]["Name"]
		vertex_prop_Seq[vertex]=G.node[node_name]["Seq"]
		vertex_prop_Cov[vertex]=G.node[node_name]["Cov"]
		vertex_prop_Color[vertex]=G.node[node_name]["Color"]
		vertex_prop_Bin[vertex]=G.node[node_name]["Bin"]
	Gt.vertex_properties["Name"] = vertex_prop_Name
	Gt.vertex_properties["Seq"] = vertex_prop_Seq
	Gt.vertex_properties["Cov"] = vertex_prop_Cov
	Gt.vertex_properties["Color"] = vertex_prop_Color
	Gt.vertex_properties["Bin"] = vertex_prop_Bin
	Dico_index_name={index:Gt.vertex_properties["Name"][Gt.vertex(index)] for index in range(Gt.num_vertices())}
	Dico_name_index={name:index for index,name in Dico_index_name.items()}
	for edges in G.edges() :
		edge1,edge2=Dico_name_index[edges[0]],Dico_name_index[edges[1]]
		Gt.add_edge(Gt.vertex(edge1),Gt.vertex(edge2))


def Save_Graph(G,Gt,file_out) :
	Handle=open(file_out+'.Nx_Pickle','w')
	pickle.dump(G, Handle)
	Handle.close()
	Gt.save(file_out+'.gt')

def Draw_Graph(Draw,Gt):
	gt.draw.graph_draw(Gt, vertex_fill_color=Gt.vertex_properties['Color'],vertex_text=Gt.vertex_properties['Name'],vertex_font_size=6,output_size=(1200,1200),output = Draw)


# contig_overlap="119"
# File_in='/home/ubuntu/DesmanExample/Example/Assembly/intermediate_contigs/k119.contigs.fa'
# contig_assignment="/home/ubuntu/DesmanExample/Example/AssignGenome/clustering_gt1000_smap_nocut.csv"
# gfa_file='/home/ubuntu/DesmanExample/Example/Assembly/intermediate_contigs/k119.contigs.gfa'
# gfa_file_raw='/home/ubuntu/DesmanExample/Example/Assembly/intermediate_contigs/k119.gfa'

contig_overlap="141"
File_in='LactoParaPhage_contig.fa'
contig_assignment=["TabC","LactoAssign.csv"]
Draw=""
# gfa_file='/home/ubuntu/DesmanExample/Example/Assembly/intermediate_contigs/k119.contigs.gfa'
# gfa_file_raw='/home/ubuntu/DesmanExample/Example/Assembly/intermediate_contigs/k119.gfa'



def main(contig_overlap,File_in,File_out,contig_assignment,Draw) :
	""" Use Megahit to generate the graph of assembly in fastg format, then use bandage to generate a gfa file, which is a compressed version. Graph data is then parsed and stored in a networkx of graph-tool datastructure. If information is availlable regarding the micro-organism or the bin, a contig come from it is coded as a color for Bandage representation and as a metadata in Networkx graph structure. """
	Dico_Contig_MO=Get_Contig_Assignment(contig_assignment)
	if File_in.split('.')[-1]=='fa' :
		fastg_file=File_in.replace('.fa','.fastg')
		gfa_file=File_in.replace('.fa','')
		os.system("megahit_toolkit contig2fastg "+contig_overlap+" "+File_in+" > "+fastg_file )
		os.system(" ".join(["Bandage reduce",fastg_file,gfa_file]))
		os.system("rm " +fastg_file)
		gfa_file=gfa_file+".gfa"
		# Add contig name and also color if contig_assignment is not empty 	
		Rewrite_gfa(Dico_Contig_MO,File_in,gfa_file)
	if File_in.split('.')[-1]=='fastg' :
		fastg_file=File_in
		gfa_file=File_in.replace('.fastg','')
		os.system(" ".join(["Bandage reduce",fastg_file,gfa_file]))
		gfa_file=gfa_file+".gfa"
	if File_in.split('.')[-1]=='gfa' :
		gfa_file=File_in
	# Extract all relevant information from gfa_file and build a Networxk graph 
	G=nx.Graph()
	Parse_gfa(gfa_file,G,Dico_Contig_MO)
	# delete redundant information, namely reverse-complement graph
	delete_redundant_comp(G,0)
	# translate graph in graph-tool object
	Gt=gt.Graph(directed=False)
	Translation_from_NX_to_Gt(G,Gt)
	Save_Graph(G,Gt,File_out)
	if Draw!="" :
		Draw_Graph(Draw,Gt)





if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("contig_overlap", help="Size of overlap between k_mer used in assembly, information needed to generate assembly graph") 
	parser.add_argument("File_in", help="Depending on task to be done can be a fasta file, a fastg or a gfa")
	parser.add_argument("Output", help="file name for graph file output")
	parser.add_argument("-C",help="contig species assignment",default="")
	parser.add_argument("-TabC",help="contig species assignment, table format",default="")
	parser.add_argument("-Draw",help="Use graph-tool to draw the assembly graph, if contig species assignment was filled, the graph is colored",default="")
	#parser.add_argument("Graph_library",help="Choose between Networkx and Graph-tool as a graph-library")
	args = parser.parse_args()
	Contig_assignment=["None",""]
	Draw=""
	if args.C :
		Contig_assignment=["C",args.C]
	if args.TabC :
		Contig_assignment=["TabC",args.TabC]
	if args.Draw :
		Draw=args.Draw
	main(args.contig_overlap,args.File_in,args.Output,Contig_assignment,Draw)

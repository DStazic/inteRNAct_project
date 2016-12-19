

###--------------------------------------------------------------
def makeGeneEntryDictionary(gene_features):
    '''take gene entry file from ncbi and make dictionary with gene names as keys and gene sequence as corresponding value
        (provide full path) '''
    
    gene_entries=open(gene_features, 'r')
    geneID = []
    geneID_seq = {}
    lines_seq=sum(1 for line in open(gene_features))
    for line in range( lines_seq):
        line_read = gene_entries.readline()[:-1]
        if ">" in line_read:
            
            geneID.append(line_read)
        else:
            if geneID[-1] not in geneID_seq.keys():
                geneID_seq[geneID[-1]]=str(line_read)
            else:
                geneID_seq[geneID[-1]]=geneID_seq[geneID[-1]]+str(line_read)
return geneID_seq
###--------------------------------------------------------------



def OneLineGenomeSeq(genome_file):
    '''take genome file and concatanate all sequence lines into one string
        (provide full path); required to match substrings that begin in one line and extend to next line'''
    genome_fasta=open(genome_file, "r")
    genome_1line=""
    for line in genome_fasta:
        if ">" not in line:
            DH10b_genome_1line=genome_1line+str(line[:-1])
    return DH10b_genome_1line

K12_genes=makeGeneEntryDictionary()
genome=OneLineGenomeSeq()
K12_sRNA = ["arcZ","chiX","cyaR", "dsrA","fnrS","gcvB", "glmZ","micA","micC", "micF","omrA","omrB","oxyS","rprA","ryhB","sgrS", "spf"]

def sRNAfind(K12_sRNA,K12_genes,genome,counter=0,all_sRNAs=[]):
    """
        Try to annotate selected K12 sRNAs for DH10b; Selection based on known sRNA:mRNA interactions characterized by Richter et al.
        1.    if sRNA name in K12_gene entry file (dict_transform_K12) take the corresponding sequence as reference for homology search in DH10b
        2.    take the first 10nt substring and search for all hits in DH10b (use regular expression module re); store coordinates for all matches in first list (list_all_coordinates_first)
        3.    Iterate over the remaining seq-10 sRNA seq and search for all hits in DH10b; ; store coordinates for all matches in second list (list_all_coordinates_second)
        4.    if a start coordinate for any hit of the next 10nt substring is within a window of 15nt next to any hit  of the first 10nt substring, store both coordinates (start/end each) in third list (list_all_coordinates_third)
        5.    For all subsequent iterations compare the next 10nt substring coordinates with the coordinates of the most recent hit coordinates in third list
        6.    After iterating over all 10nt substrings for the given K12 sRNA sequence take and store the corresponding seq in DH10b between the first 10nt substring hit START coordinate and last 10nt substring hit END coordinate
        """
    
    
    import re
    
    
    
    #make rev compl. strand  of input genome
    def RevComplGenome(genome):
        DH10b_genome_rev=""
        for base in (e for e in reversed(genome)):
            if base=="A":
                base="T"
            elif base=="T":
                base="A"
            elif base=="G":
                base="C"
            elif base=="C":
                base="G"
            
            DH10b_genome_rev+=base
        return DH10b_genome_rev
    
    
    
    
    list_all_coordinates_first=[]
    list_all_coordinates_second=[]
    list_all_coordinates_third=[]
    list_sRNAs_found=[]
    
    
    
    
    #make correct copy of K12_sRNA list; if you just do K12_sRNA_copy=K12_sRNA it will point to the same object and each time one pointer changes the other changes
    K12_sRNA_copy=list(K12_sRNA)
    
    
    """
        outer loop as generator
        """
    for sRNA_name in (e for e in K12_sRNA):
        for gene_entry in K12_genes.keys():
            try:
                
                if sRNA_name in gene_entry:
                    sRNA_seq_reference=K12_genes[gene_entry]
                    
                    #find all matches for the first 10nt substring and append to list_all_coordinates_first
                    all_matches=re.finditer(sRNA_seq_reference[:10], genome)
                    [list_all_coordinates_first.append(e.span())for e in all_matches]
                    #use list_all_coordinates_first as reference to check if the next 10nt substring is within 1-15nt range of any of the first fragment matches
                    sRNA_seq_reference=sRNA_seq_reference[10:]
                    while len(sRNA_seq_reference)>18:
                        
                        all_matches=re.finditer(sRNA_seq_reference[:10], genome)
                        [list_all_coordinates_second.append(e.span())for e in all_matches]
                        
                        # compare coordinates of second 10nt substring hits with first 10nt substring hits
                        # if any of the second substring hits coordinates is within the 15nt range of any of the first substring hits, then store both
                        
                        if list_all_coordinates_third==[]:
                            """
                                outer loop as generator
                                """
                            for frag in (e for e in list_all_coordinates_first):
                                for next_frag in list_all_coordinates_second:
                                    if next_frag[0]-frag[1]<=15 and next_frag[0]-frag[1] >= 0:
                                        
                                        
                                        list_all_coordinates_third.append(frag)
                                        list_all_coordinates_third.append(next_frag)
                    
                        # compare coordinates of subsequent 10nt substring hits with coordinates of previous 10nt substring (last item in list_all_coordinates_third)
                        # if any of the hits is within the 15nt range, then store and take this substring coordinates as reference to check against in the next iteration step
                        else:
                            for next_frag in list_all_coordinates_second:
                                if next_frag[0]-list_all_coordinates_third[-1][1]<=15 and next_frag[0]-list_all_coordinates_third[-1][1] >= 0:
                                    list_all_coordinates_third.append(next_frag)
                        #remove all coordinates from the current 10nt substring from list_all_coordinates_second; if no matches list_all_coordinates_second remains empty
                        # ==> removing all coordinates makes sure that only the coordinates of the next 10nt substring are compared with the coordinates of the preceeding 10nt substring
                        list_all_coordinates_second=[]
                        sRNA_seq_reference=sRNA_seq_reference[10:]
                
                    if len(sRNA_seq_reference)<=18:
                        all_matches=re.finditer(sRNA_seq_reference, genome)
                        [list_all_coordinates_second.append(e.span())for e in all_matches]
                        for next_frag in list_all_coordinates_second:
                            if next_frag[0]-list_all_coordinates_third[-1][1]<=15 and next_frag[0]-list_all_coordinates_third[-1][1] >= 0:
                                list_all_coordinates_third.append(next_frag)
                                    
                                    #After complete iteration of all 10nt substrings, store the start coordinate of the first hit and the end coordinate of the last hit in list_all_coordinates_third
                                    #store the DH10b seq with the cooresponding coordinates
                                    #store both as entry in list_sRNAs_found
                                    # ==> homologous K12 sRNA in DH10b
                                    #if no homologous seq found in DH10b list_all_coordinates_third remains empty; see IndexError call
                    list_sRNAs_found.append([sRNA_name,list_all_coordinates_third[0][0]+1, list_all_coordinates_third[-1][1]+1, genome[list_all_coordinates_third[0][0]:list_all_coordinates_third[-1][1]]])
                                        
                    K12_sRNA_copy.remove(sRNA_name)
                                            
                    #reset all lists for the next K12 sRNA homology search
                    list_all_coordinates_first=[]
                    list_all_coordinates_second=[]
                    list_all_coordinates_third=[]
                                                        
                    break
            #if IndexError: list index out of range for sRNA_name,list_all_coordinates_third then sRNA seq must be on rev strand
            #repeat search with DH10b genome rev strand'
            #return coordinates on complementary strand
            except IndexError:
                print ("no hit")

#Change coordinates for the sRNAs on the DH10b rev strand to match the reverse complement;
    """
    start:  surrogate for start coordinate
    end:    surrogate for end coordinate
    direct exchange --list_sRNAs_found[e][1]=len(genome)-(list_sRNAs_found[e][2]-1)-- without surogates not possible, because this will "erase" the other coordinate
    """
    if counter > 0:
        start=""
        end=""
        for e in range(len(list_sRNAs_found)):
            list_sRNAs_found[e][1]=len(genome)-(list_sRNAs_found[e][1])
            list_sRNAs_found[e][2]=len(genome)-(list_sRNAs_found[e][2]-1)
            
            start=list_sRNAs_found[e][2]
            end=list_sRNAs_found[e][1]
            
            list_sRNAs_found[e][1]=start
            list_sRNAs_found[e][2]=end
            
            all_sRNAs.append(list_sRNAs_found[e])

    else:
        [all_sRNAs.append(e) for e in list_sRNAs_found]
    
    #counter for the number of sRNAfind() function being executed
    """
        counter = 0:    first iteration over DH10b fw
        counter = 1:    second iteration over DH10b rev
        """
    counter +=1
    
    #recursive call on the remaining sRNAs in K12_sRNA_copy to be checked for homology on the DH10b reverse strand   
    """
        K12_sRNA_copy == []:                            all remaining sRNAs have been annotated on the DH10b rev strand
        counter > 1 and len(K12_sRNA_copy) != 0:        some sRNAs could not be mapped to the DH10b rev strand either due to lack of homology     
        """    
    if K12_sRNA_copy==[] or (counter > 1 and len(K12_sRNA_copy) != 0): 
        return all_sRNAs           
    
    
    
    revStrandSearch=sRNAfind(K12_sRNA_copy,K12_genes, RevComplGenome(genome),counter,all_sRNAs)
    return revStrandSearch

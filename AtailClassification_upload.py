
def AtailClassification(Path_fastQ_InQuotes,output_name_InQuotes,Atail_length=5):
    '''Classify reads into potential interaction sites versus no interaction site (no proximity ligation)
       Potential interaction sites contain internal 8-20nt A-tails 
       
       potential interaction site = NNNNAAAANNNN (N=any rRNA base) or NNNNAAAANNNNAAAA
       no interaction site = NNNN or NNNNAAAA
       
       Path_fastQ_InQuotes = full path of fastQ file (use quotes for path in function)
       output_name_InQuotes = name of output file (use quotes for file name in function)
       Atail_length = kmer size; 5nt default
       
       
       1. take fastQ file as input
       2. extract the sequence for each read
       3. define kmer size (nt length) and calculate A (Adenosine) percentage in kmer
       4. start at 3'end of read and move kmer along the sequence by one nucleotide
       5. for every subsequent position calculate A percentage in kmer
       5.1  try to find following pattern = kmer: A% < 100 kmer: A% =100 kmer: A% < 100
            --> potential interaction site
       5.2  else no interaction site
        '''
    
    # generator to count the number of lines in input fastQ file
    # generate a list of lines to be skipped by pandas
    # thus, read only sequence lines from fastQ file
    '---> every 4th line starting at 2; In python starting with line 1, because line 0 = first line!!!'
    num_lines = sum(1 for l in open(Path_fastQ_InQuotes))
    idx_all_lines=range(0,num_lines)
    idx_skip_lines=[x for x in idx_all_lines if (x+3)%4 !=0]
    
    import pandas
    
    fastQ=pandas.read_csv(Path_fastQ_InQuotes, header=None, skiprows=idx_skip_lines)
    output_name=open("/data/inteRNAct/raw/20160829/Term.Transf.Efficiency/"+output_name_InQuotes,"w")
    for e in fastQ[0]:
        ### check first kmer for ------------------------------
        
        first_kmer=e[::-1][0:Atail_length]
        count_A_first_kmer=0
        for base_first_kmer in first_kmer:
            if base_first_kmer=="A":
                count_A_first_kmer+=1
        percent_A_first_kmer=(count_A_first_kmer/float(Atail_length))*100
        percent_next_kmer=None
        
        ### no proximity ligation -------------------------------------------------------------------------------------------
        if percent_A_first_kmer==100.0:
            count_A=0
            start_Atail=None
            start_internal_Atail=None
            end_internal_Atail=None
            
            for iteration in range(1,len(e[:-Atail_length])):                    #-5 iterations, because last len(kmer) are not iterated over; len(e)-len(kmer) iterations left
                next_kmer=e[::-1][iteration:Atail_length+iteration]          #next 5-kmer starting from base 2
                for base in next_kmer:
                    if base=="A":
                        count_A+=1
                percent_next_kmer=(count_A/float(Atail_length))*100
                count_A=0
                
                if percent_next_kmer<100.0:
                    if start_Atail==None:                                          #calculate the end of A-tail only for the first consecutive kmer <100%!!
                        start_Atail=(len(e)-iteration-len(next_kmer)-1)+1
                    if end_internal_Atail != None and start_internal_Atail==None:
                        start_internal_Atail=(len(e)-iteration-len(next_kmer)-1)+1    #calculate the end of internal A-tail only for the first consecutive kmer =100%!!
                '+len(next_kmer)-1:    -1 because the last base of kmer does not belong to A-tail'
                '+1:                   +1 because start_Atail equation would otherwise return the position of the base that is upstream of A-tail'
                if percent_next_kmer==100.0 and start_Atail != None:
                    if end_internal_Atail==None:                               #calculate the start of internal A-tail only for the first consecutive kmer =100%!!
                        end_internal_Atail=len(e)-iteration
            
            if start_internal_Atail==None:         #terminal transferase class
                
                output_name.write("no proximity ligation"+"\t"+"start:"+" "+str(start_Atail)+"\t"+"end:"+" "+str(len(e))+"\t"+str(e)+"\n")
    else:
        
        output_name.write("no proximity ligation"+"\t"+"start:"+" "+str(start_internal_Atail)+"\t"+"end:"+" "+str(len(e))+"\t"+str(e)+"\n")
        ### potential interaction site -----------------------------------------------------------------------------
        elif percent_A_first_kmer<100.0:
            count_A=0
            start_Atail=None
            end_Atail=None
            for iteration in range(1,len(e[:-Atail_length])):
                next_kmer=e[::-1][iteration:Atail_length+iteration]
                for base in next_kmer:
                    if base=="A":
                        count_A+=1
                percent_next_kmer=(count_A/float(Atail_length))*100
                count_A=0
                
                if percent_next_kmer==100.0:
                    # similar to no proximity ligation
                    # keep first kmer position of kmer=100%A --> this is the start position of the polyA tail
                    if end_Atail==None:
                        end_Atail=len(e)-iteration
        
                elif percent_next_kmer<100.0 and end_Atail !=None:
                    start_Atail=(len(e)-len(e[:iteration+len(next_kmer)-1]))+1
                    '+len(next_kmer)-1:    -1 because the last base of kmer does not belong to A-tail'
                    '+1:                   +1 because start_Atail equation would otherwise return the position of the base that is upstream of A-tail'
                    if len(e)-end_Atail>=10:
                        'filter reads with at least 10nt after A-tail (treshhold for min RNA length)'
                        output_name.write("Potential interaction site"+"\t"+"start:"+" "+str(start_Atail)+"\t"+"end:"+" "+str(end_Atail)+"\t"+str(e)+"\n")
                        break
                else:
                break
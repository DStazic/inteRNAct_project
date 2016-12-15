#start at 3' end
#if first kmer=100% ----> 3'A-tail

#2. continue until next(second) kmer <100% ----->end(actually start) of 3'A-tail
#3. continue until next(third) kmer=100% ----->start(actually end) of term.transf A-tail
#4. continue until next kmer(fourth)<100% ----->end(actually start) of term.transf A-tail
#5. if end of read is reached and third and fourth kmer=None-----> read belongs to 3'Atail class
#   if end of read is reached and third and fourth kmer not None-----> read belongs to terminator transferase class

#if first kmer<100% ----> potentially terminator transferase read

#2. continue until next(second) kmer =100% ----->start(actually end) of term.transf A-tail
#3. continue until next kmer(third)<100% ----->end(actually start) of term.transf A-tail
#4. if end of read is reached and second and third kmer=None-------> read belongs to standard class
#   if end of read is reached and second and third kmer not None-------> read belongs to terminator transferase class






def AtailClassification(Path_fastQ_InQuotes,output_name_InQuotes,Atail_length=5,sep_="\t"):
    # generator to count the number of lines in input fastQ file
    # generate a list of lines to be skipped by pandas
    # thus, read only sequence lines from fastQ file
    '---> every 4th line starting at 2; In python starting with line 1, because line 0 = first line!!!'
    num_lines = sum(1 for l in open(Path_fastQ_InQuotes))
    idx_all_lines=range(0,num_lines)
    idx_skip_lines=[x for x in idx_all_lines if (x+3)%4 !=0]
    
    import pandas
    Path_fastQ_InQuotes=str(Path_fastQ_InQuotes)
    fastQ=pandas.read_csv(Path_fastQ_InQuotes, header=None,sep="\t", skiprows=idx_skip_lines)
    output_name=open("/data/inteRNAct/raw/20160829/Term.Transf.Efficiency/"+str(output_name_InQuotes),"w")
    for e in fastQ[0]:
        ### for Adapter trimmed reads check if first kmer is potentially solexa A-tail------------------------------
        ### best take smaller kmer -5nt-for this initial step
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
                    # similar to 3'A-tail path
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
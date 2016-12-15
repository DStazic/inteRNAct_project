
######------------------------------------------------------------------------------------------------------
# load reference table with gene entries (column 3) and corresponding start- (column 1) and end coordinates (column 2)

reference =read.table("/Users/damirvana/Arbeit_Uni/Postdoc_AG_Voß/Sequenzierung/E.coli_DH10B_K12_features_gene.coordinates.NAMES.txt")


# load table with split reads (mapped as split reads using tophat or segemehl)
# 
split_reads = read.table("/Users/damirvana/Arbeit_Uni/Postdoc_AG_Voß/Sequenzierung/kamo_251/splicesites_kamo_251.bed")



GeneIdentify=function(reference, split_reads) {
######------------------------------------------------------------------------------------------------------
# make a copy of the split_reads table for modifications  
split_reads_copy=split_reads  

######------------------------------------------------------------------------------------------------------
# split the reference list with all gene entries into sets of 100 genes (43 sets for E.coli DH10b)
# reference list was initially sorted for gene start/end coordinates in ascending order 
# calculate the min() and max() coordinate for all gene sets 
# coordinate for each mapped read will be crosschecked if it is within any of the gene set min/max boundaries
# loop over the correct gene set and reassign the mapped read start position to the corresponding gene annotation 
'-----> dividing the total set of gene entries into sets of 100 genes speeds up computation'
'-----> for each read (millions in total!!) loop over 1. (max)43 sets and over 2. (max)100 entries in one set instead of 4300 gene entries'

split_reference_list=list()   
split_reference_part=data.frame()^
for (i in 1:nrow(reference)) {
  split_reference_part=rbind(split_reference_part, reference[i,],stringsAsFactors=FALSE)
  
  if (i%%100==0) {
    name <- paste('index:',i,sep='')
    split_reference_list[[name]]=split_reference_part
    split_reference_part=data.frame()
  }
}
# calculate min() and max() coordinate for all groups; this can be used if coordinates in split reads 
# are within min()/max() bounderies of any of these groups. If TRUE then search only in this group for hit
# here, first search in approx. 40 (list with min/max for every 100 genes set), then search in approx. 100 (list with sets of 100 genes). 
# (In total 140 vs. 4300 (list with all genes) entries to check.

lower=c()
upper=c()
gene_set.boundaries=data.frame()

for (i in 1:length(split_reference_list)) {
  lower=min(split_reference_list[[i]]$V1)
  upper=max(split_reference_list[[i]]$V2)
  gene_set.boundaries=rbind(gene_set.boundaries,c(lower,upper), stringsAsFactors=FALSE)
}

colnames(gene_set.boundaries)=c("Coordinate_start","Coordinate_end")
rownames(gene_set.boundaries)=names(split_reference_list)

######------------------------------------------------------------------------------------------------------

summary_split_reads=data.frame()
coord1_ref=c()
coord1_split=c()
name1=c()

for (split in 1:nrow(split_reads_copy)) {
  
    if (split_reads$V1[split]=="gi|170079663|ref|NC_010473.1|") {
      for (ref in 1:nrow(gene_set.boundaries)) {
        if (gene_set.boundaries$Coordinate_start[ref]<=split_reads$V2[split] & gene_set.boundaries$Coordinate_end[ref]>=split_reads$V2[split]) {
          list_to_look_up_first=as.data.frame(split_reference_list[[ref]])
                                     
          for(i in 1:nrow(list_to_look_up_first)) {
            if (list_to_look_up_first$V1[i]<=split_reads$V2[split] & list_to_look_up_first$V2[i]>=split_reads$V2[split]) {
              
              #here copy list                            
              split_reads_copy$V2[split]=as.character(list_to_look_up_first$V3[i])
              list_to_look_up_first=data.frame()
              #print(i)
                                          
              break
                                          
          }
          }
          
        }
        ###
        if (gene_set.boundaries$Coordinate_start[ref]<=split_reads$V3[split] & gene_set.boundaries$Coordinate_end[ref]>=split_reads$V3[split]) {
          list_to_look_up_second=as.data.frame(split_reference_list[[ref]])
          
          for(i in 1:nrow(list_to_look_up_second)) {
            if (list_to_look_up_second$V1[i]<=split_reads$V3[split] & list_to_look_up_second$V2[i]>=split_reads$V3[split]) {
              
              #here copy
              split_reads_copy$V3[split]=as.character(list_to_look_up_second$V3[i])
              list_to_look_up_second=data.frame()
              print(i)
              
              break
              
            }
          }
          break
        }
        
    }
    
  }
}

return(split_reads_copy)

}


write.table(split_reads_copy,file="/Users/damirvana/Arbeit_Uni/Postdoc_AG_Voß/Sequenzierung/DH10b_trivial.gene.annotations_splicesites_kamo_251.txt",row.names=FALSE,sep="\t",col.names = FALSE)

write.table(split_reads_copy,file="/Users/damirvana/Arbeit_Uni/Postdoc_AG_Voß/Sequenzierung/DH10b_trivial.gene.annotations_splicesites_kamo_250.txt",row.names=FALSE,sep="\t",col.names = FALSE)


##########






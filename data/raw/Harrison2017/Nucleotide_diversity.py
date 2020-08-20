import sys


#my_vcf = open(sys.argv[1],'r')

my_vcf =open("Mlupulina_snps.vcf",'r')


print "chrom", " " , "samp1", " ", "samp2"," ", "num_sites_where_both_havedata", " ","var_sites"," ", "N_both_het"," ",   "N_hom_het"," ", "N_hom_alt_hom_ref", " ", "N_hom_hom",  " ", "pairwise_geno"


def CalcWind(lastChrom, N_both_het, N_hom_het, N_hom_alt_hom_ref,  N_hom_hom, pairwise_geno,site_num, var_num):

    for key in pairwise_geno.keys():
        #print key #keys' are sample names
        # v is the dictionary of stuff 
        for keys in pairwise_geno[key].keys():
            print lastChrom, key, keys, site_num[key][keys], var_num[key][keys],N_both_het[key][keys], N_hom_het[key][keys], N_hom_alt_hom_ref[key][keys],  N_hom_hom[key][keys], pairwise_geno[key][keys]



indIDs =['AIL_1_3_mel_reorder.bam','AV_1_1_mel_reorder.bam','AV_1_4_mel_reorder.bam','BLU1_2_mel_reorder.bam','BLU1_6_mel_reorder.bam','BUF_1_3_mel_reorder.bam','CAL_1_2_mel_reorder.bam','CAL_1_6_mel_reorder.bam','CAL_2_1_mel_reorder.bam','CAL_2_2_mel_reorder.bam','COO_1_4_mel_reorder.bam','COO_1_5_mel_reorder.bam','DE_1_2_mel_reorder.bam','DE_1_3_mel_reorder.bam','DE_1_4_mel_reorder.bam','DUN_1_2_mel_reorder.bam','DUN_1_6_mel_reorder.bam','FOR_1_4_mel_reorder.bam','GA_1_5_mel_reorder.bam','GIL_1_1_mel_reorder.bam','GIL_1_2_mel_reorder.bam','HA_1_2_mel_reorder.bam','HA_2_4_mel_reorder.bam','HL_1_4_mel_reorder.bam','HL_1_6_mel_reorder.bam','KA_1_3_mel_reorder.bam','KIT_1_6_mel_reorder.bam','KSR_1_1_mel_reorder.bam','MTC_1_3_mel_reorder.bam','NEY_1_5_mel_reorder.bam','NEY_1_7_mel_reorder.bam','NEY_1_8_mel_reorder.bam','NH_1_5_mel_reorder.bam','PAN_1_2_mel_reorder.bam','PAN_1_3_mel_reorder.bam','PAN_1_6_mel_reorder.bam','PAR_1_1_mel_reorder.bam','PA_1_2_mel_reorder.bam','PA_1_6_mel_reorder.bam','PA_1_7_mel_reorder.bam','PA_1_8_mel_reorder.bam','PT_1_10_mel_reorder.bam','PT_1_4_mel_reorder.bam','PT_1_5_mel_reorder.bam','PUS_1_3_mel_reorder.bam','PUS_1_5_mel_reorder.bam','SC_1_3_mel_reorder.bam','SC_1_4_mel_reorder.bam','SIN_1_1_mel_reorder.bam','SIN_1_2_mel_reorder.bam','TO_3_1_mel_reorder.bam','TO_3_2_mel_reorder.bam','UBR_1_2_mel_reorder.bam','UBR_1_3_mel_reorder.bam']
                    
                    
                    
def ResetDicts():
    site_num={}
    var_num={}

    pairwise_geno={}
    N_both_het={}
    N_hom_het={}
    N_hom_alt_hom_ref={}
    N_hom_hom={}

    indIDs2=indIDs
    for x in indIDs:
        #pairwise_pi is a dictionary where the keys aresample_ID names, and then each sample name is associated with a dictionary of all pairwise comparisons among samples.
        site_num[x]={}
 
        pairwise_geno[x]={}
        N_both_het[x]={}
        N_hom_het[x]={}
        N_hom_alt_hom_ref[x]={}
        N_hom_hom[x]={}
        var_num[x]={}
  
        indIDs2=[item for item in indIDs2 if item != x]

        for k in (indIDs2):

            pairwise_geno[x][k]=0.0
            site_num[x][k]=0.0
            N_both_het[x][k]=0.0
            N_hom_het[x][k]=0.0
            N_hom_alt_hom_ref[x][k]=0.0
            N_hom_hom[x][k]=0.0
            var_num[x][k]=0.0

    return N_both_het, N_hom_het, N_hom_alt_hom_ref,  N_hom_hom, pairwise_geno,  site_num,var_num


N_both_het, N_hom_het, N_hom_alt_hom_ref,  N_hom_hom, pairwise_geno, site_num,var_num = ResetDicts()    




start=-1
new_line=[]
sample2={}
last_pos=0
variant_sites=0




lastChrom="Chromosome"
count=0
for line in my_vcf:
    count+=1

    if count==100:
        CalcWind(lastChrom, N_both_het, N_hom_het, N_hom_alt_hom_ref,  N_hom_hom, pairwise_geno,  site_num,var_num)


    if line[0]=="#" and line[0:2]!="#C":
        continue

#once we get to line with sample names output pariwise comparisons of samples, 9 =position that python starts reading sample names      
    if line[0:2]=="#C":
        sline=line.split()
        indIDs = sline[9:]
        N_both_het, N_hom_het, N_hom_alt_hom_ref,  N_hom_hom, pairwise_geno,  site_num, var_num = ResetDicts()
        continue
 
    

    
    sline=line.split()
    
    #if sline[0]=="pseudo0":
    #    continue    
#_N is going to be number of samples at that site
  #  _N= sum(1 for x in sline if my_condition(x))
    currentChrom=sline[0]
    
    
    if currentChrom!=lastChrom:
        sys.stderr.write("processing %s  1 \n" %(currentChrom)) 
        
        CalcWind(lastChrom, N_both_het, N_hom_het, N_hom_alt_hom_ref,  N_hom_hom, pairwise_geno, site_num, var_num)

        N_both_het, N_hom_het, N_hom_alt_hom_ref,  N_hom_hom, pairwise_geno, site_num,var_num = ResetDicts()

#For Tia we have to say which columns are genotype data
    genos=sline[9:]    

    for k in range(0,len(indIDs)):
       # print k
        samp1=indIDs[k]
        
        #if this sample has no data then skip doing all the comparisons for it
        if genos[k][0]!=".":
            print genos[k][0]
            #this will compare the first sample to all others in the row
            for i in range(k+1,len(indIDs)):
                
                samp2=indIDs[i]

        
                if genos[i][0]!=".":
           #both samples have to have data (meaning both are not  N for this section to compute)
                    site_num[samp1][samp2]+=1
                    
#############################                    

                #both sames are homozygote for same alleles
                if (genos[k][0]=="0" and genos[i][0]=="0") or (genos[k][0]=="1" and genos[i][0]=="1"):
                    N_hom_hom[samp1][samp2]+=1
                #samples are alt homoz to each other
                elif (genos[k][0]=="0" and genos[i][0]=="1") or (genos[k][0]=="1" and genos[i][0]=="0"):
                    N_hom_alt_hom_ref[samp1][samp2]+=1
                    pairwise_geno[samp1][samp2]+=1

                    print   pairwise_geno[samp1][samp2]


    lastChrom=sline[0] 





CalcWind(lastChrom, N_both_het, N_hom_het, N_hom_alt_hom_ref,  N_hom_hom, pairwise_geno, site_num, var_num)

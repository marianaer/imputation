#!/usr/bin/env python


import datetime
import string
from sys import argv
import numpy as np
import  re
import math


startTime = datetime.datetime.now()

np.set_printoptions(threshold=np.inf)


# How to use: we can chose which populations to convert to vcf (Pop1,2,3 or all of them). Coverage and maf values can be set to the user's preference
# python fsc2vcf_lowCovs.py  arpfile.arp outfilename vcf pop1 pop2 pop3 cov=1 maf=0.05


#############################
### FUNCTION TO SIMULATE MISSING RATE ###
#############################


def phred2prob(x):
    return 10.0**(-x/10.0) #error

def prob2phred(x):
    return -10*math.log10(x)

all_dic={}
all_dic['A']=0
all_dic['C']=1
all_dic['G']=2
all_dic['T']=3

def geno_caller_3GT(X,ref,alt,all_dic):
    #diploid caller assuming that only assesses likelihood for three possible genotypes (ref/ref,ref/alt,alt/alt)           
    GL=[0.0,0.0,0.0]

    count=0
    for g in range(len(X)):
        if all_dic.has_key(X[g][0])==False:
            continue
        err=phred2prob(X[g][1]) #el input es el error en phred
        tru=1-phred2prob(X[g][1])
      
        if X[g][0]==ref:
            GL[0]=GL[0]+math.log10(tru)
            GL[1]=GL[1]+math.log10((tru+err)/2)
            GL[2]=GL[2]+math.log10(err)           
        elif X[g][0]==alt:
            GL[0]=GL[0]+math.log10(err)
            GL[1]=GL[1]+math.log10((tru+err)/2)
            GL[2]=GL[2]+math.log10(tru)  
        else:
            GL[0]=GL[0]+math.log10(err) 
            GL[1]=GL[1]+math.log10(err)
            GL[2]=GL[2]+math.log10(err)
        count+=1

    if count==0:
        GL=[-9.0,-9.0,-9.0]
    return GL 




filein = argv[1]
filenameout= argv[2]

#fileout = open(filenameout+'.ped','w')

loci=[]
ind_counter=0
ind_list=[]
chrs_sites=[]
chrs_pos=[]
chr_ord=[]
sample_s=[]

#To store the positions of each snp of each chromosome
with open(filein,'r') as file: 
    for line in file:
        if 'SampleSize' in line:
            sample_s.append(int((line.split('='))[1]))
        if '#Number of independent chromosomes:' in line:
            n_chrs=int(line.split(':')[1])
        if 'polymorphic positions on chromosome' in line:
            line=line.replace('#', '')
            ch=int((line.split('polymorphic')[1]).split('chromosome')[1]) #prints the number of each chr
            print(ch)
            line=next(file)
            line=line.replace('#','')
            line=line.replace(' ', '')
            line=line.replace("'", "")
            line=line.split(",")
            print(len(line))
            chr_ord.append(np.repeat(ch,len(line))) #how many sites each chr has
            line=[int(i) for i in line]
            chrs_pos.append(line) #appends all the positions of the sites of that chr (a list of lists)
chrs_pos = np.array([item for sublist in chrs_pos for item in sublist]) # converts the list of lists into one list. stores all the positions of the sites
chr_ord= np.array([item for sublist in chr_ord for item in sublist]) 

n_chrs=np.unique(chr_ord)
n_snps=len(chrs_pos)

print('Total number of polymporphic sites:', n_snps)



len_pop1_hap=sample_s[0]
len_pop2_hap=sample_s[1]
len_pop3_hap=sample_s[2]




### PARSING ARP FILE

file = open(filein,'r')

x = file.readline()

allele_dic = {0:'A', 1:'T', 2:'C', 3:'G'}


while x!='': # while theres a line
    x=x.replace('\r','') # remove
    if '[Profile]' in x[:-1]: # if its the profile line
        x = file.readline() #  read next two lines
        x = file.readline() # 
        x=x.replace('\r','') 
        nb_samps=int(string.split(x[:-1],'=')[1]) # read how many samps

        sample_nb=0 # counts # of samples (populations)

################################################################################

        
    if 'Reporting status of a maximum of' in x[:-1]:
        loci=[]
        n_snps=int(string.split(x[:-1],)[6])#se queda solo con el # y lo guarda
        loci.append(n_snps)
	chrs_pos=chrs_pos[:n_snps] #keep only as much positions as sites we have
	print('n_sps: ', n_snps) 


################################################################################

    if 'SampleName=' in x[:-1]: 
        x = x.replace('\r', '') 
        x = x.replace(' ','')
        t22 = x.rsplit('=', -1)
        t23 = str(t22[1:])
        t23 = t23.replace('\\n', '')
        #t23.translate(str.maketrans('', '', string.punctuation))
        t23 = t23.translate(None, string.punctuation) # name of the sample

    
	
# Keep # of inds
    if 'SampleSize=' in x[:-1]:
    	x = x.replace('\r', '')
    	inds = (string.split(x[:],"="))
    	inds = int(inds[1])
    	nb_ind=sum(sample_s)

        if ind_counter==0:
            geno_array=np.zeros((n_snps,nb_ind),dtype='int32')
################################################################################


    if 'SampleData= {' in x[:-1] : 
        sample_nb+=1 # Adds to sample counter
        #print(x)
        x=file.readline() # reads next line
        x=x.replace('\r','') 


        while x!='\n' : # as long as its reading the line
           

        ## Parse genotype line


            x=x.split()
            x[2]=list(x[2])
            ini=x[0]+'\t'+x[1]+'\t '
            x=ini+str(list(x[2]))
            x=x.replace(",","") 
            x=x.replace("[","")
            x=x.replace("]","")
            x=x.replace("\'","")
            x=x.replace('"', '')
            x=x.replace('\r','')

            y1=string.split(x)
            x=file.readline() 
            x=x.split()
            x[2]=list(x[2])
            ini=x[0]+'\t'+x[1]+'\t '
            x=ini+str(list(x[2]))
            x=x.replace(",","") 
            x=x.replace("[","")
            x=x.replace("]","")
            x=x.replace("\'","")
            x=x.replace('"', '')
            x=x.replace('\r','') 
            y2=string.split(x) 
            ttt25 = y1[0] #Guarda ID 1_1 etc


            if len(y1)==len(y2): 
                genos1=y1[2:] # stores genotype data (los 0 y 1)
                genos1=genos1[0:n_snps]
                genos2=y2[2:] # stores genotype data (los 0 y 1)
                genos2=genos2[0:n_snps]		

            elif len(y1)==len(y2)+2: 
                genos1=y1[2:] #Guarda solo datos de geno (los 0 y 1)
                genos2=y2[:] #Guarda toda la linea


                    
            geno_array[:,ind_counter]=np.array(genos1) # stores genotype data
            ind_counter+=1 
         
            geno_array[:,ind_counter]=np.array(genos2)
            ind_counter+=1
            #print(ind_counter)


            if x=='\n': # If its not sample data anymore
                break # stop looping


            x=file.readline()
            x=x.replace('\r','')

    x=file.readline()
#fileout.close() 


mask=np.zeros(n_snps,dtype='int32') # create mask array for triallelic snps


for g in range(n_snps):
    if len(np.unique(geno_array[g]))!=2: # if the site doesent have only two alleles
        mask[g]=1

usable_snps=np.where(mask==0) # Keep index of only biallelic snps
geno_usable=geno_array[usable_snps] #array with only biallelic snps
chr_ord=chr_ord[usable_snps] #array with chromosomes of only biallelic snps
chrs_pos=chrs_pos[usable_snps] #positions of biallelic snps

print('After removing triallelic sites:',len(chr_ord))


#####################
#### MAF filter #####
#####################

maf_arg = re.compile("maf*") 
#look for "missing" in argv array
maf = filter(maf_arg.match, argv)
maf=str(maf)
maf=maf.split("=")
maf=str(maf[1])
maf=maf.replace("'", '')
maf=maf.replace("]", '')
maf=maf #Stores missing %

#nb_ind = # of alleles
counts_list=[] # keeps the index of sites

for i in range(geno_usable.shape[0]): #goes through each snp

    c=((np.unique(geno_usable[i,:len_pop1_hap], return_counts=True))[1]) # how many counts we have of each allele . We should filter by MAF on the ref Pop (pop1)

    if min(c) > (float(maf)*len_pop1_hap): #here it should be the MAF * the size of Pop1 *2 (bc 2 copies) so haploid length works
        counts_list.append(i) #add the index of this site to the list

geno_usable=geno_usable[counts_list] #indices<<<< of usable sites according to MAF

g1=np.vectorize(allele_dic.get)(geno_usable[:,::2]) ## Odd indexes (1,3,5...) genotypes form ind 1_1, 1_3...

g2=np.vectorize(allele_dic.get)(geno_usable[:,1::2]) ## Even indexes (2,4,6...) genotypes from 1_2, 1_4...


# store the genomic positions filtered by maf
usable_pos=chrs_pos[counts_list] # Usable positions filtered by maf
chr_ord=chr_ord[counts_list] #chromosome of each site

print('AFter MAF filtering:', len(chr_ord))


 #####################
 ### VCF #############
 #####################

inds=[]
for s, l  in enumerate(sample_s[0:3]): #we only care about pop1 and pop2
    for i in range(0, l):
        if i % 2 == 0:
            inds.append(str(s+1)+'_'+str(i+1)+'_'+str(s+1)+'_'+str(i+1))


inds_pop1=inds[:len_pop1_hap/2]
inds_pop2=inds[len_pop1_hap/2:(len_pop1_hap/2+len_pop2_hap/2)]
inds_pop3=inds[(len_pop1_hap/2+len_pop2_hap/2):]

inds_p2=inds_pop2

phase_arr=np.zeros((len(geno_usable), (len_pop1_hap+len_pop2_hap+len_pop3_hap)/2), dtype='object') # Divide by 2 to output pop1 and pop2 as diploid organisms (we have 3 pops)


for i in range((len_pop1_hap+len_pop2_hap+len_pop3_hap)/2): # for each individual in Pop1 and Pop2
    for g in range(len(geno_usable)): # for each snp

        g1_allele= g1[g,i]        
        g2_allele= g2[g,i]
        ref_alt=np.unique(np.concatenate((g1[g,:], g2[g,:])))
        if len(ref_alt)<2: #if the site is monomorphic. 
            phase_arr[g,i]='0|0'

        else:        
            if g1_allele==ref_alt[0] and g2_allele==ref_alt[0]: # if they're the same BUT the site is  not monomorphic= hom ref
                phase_arr[g,i]='0|0'


            elif g1_allele==ref_alt[1] and g2_allele==ref_alt[1]: # if they're the same BUT the site is  not monomorphic= hom alt
                phase_arr[g,i]='1|1'

            elif g1_allele==ref_alt[0] and g2_allele==ref_alt[1]: # heterozygous site
                phase_arr[g,i]='0|1'

            elif g1_allele==ref_alt[1] and g2_allele==ref_alt[0]: # heterozygous site
                phase_arr[g,i]='1|0'


# CONVERT TO VCF

if 'vcf' in argv:

	d = datetime.datetime.now()
	date=(d.strftime("%Y"),d.strftime("%x").split('/')[0],d.strftime("%x").split('/')[1])
	date=str(date)
	date=date.replace("'", '')
	date=date.replace(",", '')
	date=date.replace("(", '')
	date=date.replace(")", '')
	date=date.replace(" ", '')
	date=('##fileDate='+date)


p = re.compile("pop2") #look for "pop2" in argv array. this means we want phased true pop2
p2_arg = filter(p.match, argv)

pp = re.compile("pop1") # if we want to output the reference pop
p1_arg=filter(pp.match, argv)

ppp = re.compile("pop3") # If we want to ouotput the outgroup
p3_arg=filter(ppp.match, argv)

if p1_arg or p2_arg or p3_arg:

    for i in range(1,4):
        for j in n_chrs :
            with open(filenameout+'_'+maf.split(".")[1]+'maf'+'_Pop'+str(i)+'_chr'+str(j)+'_phased'+'.vcf', 'w') as chr_file:
                chr_file.write('##fileformat=VCFv4.1\n')
                chr_file.write(date+'\n')
                chr_file.write('##contig=<ID='+str(j)+'>\n')
                chr_file.write('##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n')
                chr_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                header=('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')

                inds_pop1=str(inds_pop1)
                inds_pop1=inds_pop1.replace("'", '')
                inds_pop1=inds_pop1.replace("[", '')
                inds_pop1=inds_pop1.replace("]", '')
                inds_pop1=inds_pop1.replace(",", '')
                inds_pop1=inds_pop1.replace(" ", '\t')

                inds_pop2=str(inds_pop2)
                inds_pop2=inds_pop2.replace("'", '')
                inds_pop2=inds_pop2.replace("[", '')
                inds_pop2=inds_pop2.replace("]", '')
                inds_pop2=inds_pop2.replace(",", '')
                inds_pop2=inds_pop2.replace(" ", '\t')


                inds_pop3=str(inds_pop3)
                inds_pop3=inds_pop3.replace("'", '')
                inds_pop3=inds_pop3.replace("[", '')
                inds_pop3=inds_pop3.replace("]", '')
                inds_pop3=inds_pop3.replace(",", '')
                inds_pop3=inds_pop3.replace(" ", '\t')


                chr_sites=np.where(chr_ord==j) # Keeps the sites in this chr
                geno_pos=usable_pos[chr_sites]
                phase_chr=phase_arr[(chr_sites)] # the phasing of this site for all individuals
                chr_genos=g1[chr_sites] 


                if i == 1 and p1_arg: # Pop1
                    print 'Generating Pop1...'
                    chr_file.write(header+inds_pop1+'\n')# only outputs individuals from pop1


                    for l in range(len(geno_pos)): # For each snp in this chr
                        pos=str(geno_pos[l]) #keep theposition
                        ref=str(np.unique(chr_genos[l,:])[0]) #keep the ref allele

                        if len(np.unique(chr_genos[l,:]))==2: 
                            alt=str(np.unique(chr_genos[l,:])[1]) #Keep the alt allele

                            out= str(j)+ '\t'+pos+'\t'+'snp_'+str(j)+pos+'\t'+ref+'\t'+alt+'\t.\t.\tPR\t'+'GT'+'\t'+str(str(phase_chr[l, :len_pop1_hap/2]).split())+'\n'
                            out=out.replace('[','')
                            out=out.replace(']','')
                            out=out.replace('"','')
                            out=out.replace("'",'')
                            out=out.replace(",",'')
                            out=out.replace(" ",'\t')
                            chr_file.write(out)

                elif i == 2 and p2_arg: # Pop2    
                    print('Generating Pop2...')        
                    chr_file.write(header+inds_pop2+'\n')# only outputs individuals from pop2


                    for l in range(len(geno_pos)): # For each snp
                        pos=str(geno_pos[l]) #keep theposition
                        ref=str(np.unique(chr_genos[l,:])[0]) #keep the ref allele

                        if len(np.unique(chr_genos[l,:]))==2: 
                            alt=str(np.unique(chr_genos[l,:])[1]) #Keep the alt allele

                            out= str(j)+ '\t'+pos+'\t'+'snp_'+str(j)+pos+'\t'+ref+'\t'+alt+'\t.\t.\tPR\t'+'GT'+'\t'+str(str(phase_chr[l, len_pop1_hap/2:(len_pop1_hap/2+len_pop2_hap/2)]).split())+'\n'
                            out=out.replace('[','')
                            out=out.replace(']','')
                            out=out.replace('"','')
                            out=out.replace("'",'')
                            out=out.replace(",",'')
                            out=out.replace(" ",'\t')
                            chr_file.write(out)


                elif i == 3 and p3_arg: # Pop3    
                    print('Generating Pop3 (outgroup)...')        
                    chr_file.write(header+inds_pop3+'\n')# only outputs individuals from pop3


                    for l in range(len(geno_pos)): # For each snp
                        pos=str(geno_pos[l]) #keep theposition
                        ref=str(np.unique(chr_genos[l,:])[0]) #keep the ref allele

                        if len(np.unique(chr_genos[l,:]))==2: 
                            alt=str(np.unique(chr_genos[l,:])[1]) #Keep the alt allele

                            out= str(j)+ '\t'+pos+'\t'+'snp_'+str(j)+pos+'\t'+ref+'\t'+alt+'\t.\t.\tPR\t'+'GT'+'\t'+str(str(phase_chr[l, (len_pop1_hap/2+len_pop2_hap/2):]).split())+'\n'
                            out=out.replace('[','')
                            out=out.replace(']','')
                            out=out.replace('"','')
                            out=out.replace("'",'')
                            out=out.replace(",",'')
                            out=out.replace(" ",'\t')
                            chr_file.write(out)
                          


### If user specified to simulate coverage:
r = re.compile("cov*") #look for "missing" in argv array
m = filter(r.match, argv)
if m: #If user did specify a coverage:
    print('Generating ancient individuals...')
    m=str(m)
    m=m.split("=")
    m=str(m[1])
    m=m.replace("'", '')
    m=m.replace("]", '')
    m=m #Stores coverage %

    d = datetime.datetime.now()
    date=(d.strftime("%Y"),d.strftime("%x").split('/')[0],d.strftime("%x").split('/')[1])
    date=str(date)
    date=date.replace("'", '')
    date=date.replace(",", '')
    date=date.replace("(", '')
    date=date.replace(")", '')
    date=date.replace(" ", '')
    date=('##fileDate='+date)


    
    for indx, kk in enumerate(inds): 

        if kk in inds_p2:

            i1=np.where(np.array(inds)==inds_p2[0]) #keeps the index of the first ind of pop2
            i1=i1[0][0]


            for j in n_chrs:
                with open(filenameout+'_'+maf.split(".")[1]+'maf_'+str(m)+'cov_'+'ind_'+kk+'_chr'+str(j)+'_unphased'+'.vcf', 'w') as chr2_file:
                    
                    chr2_file.write('##fileformat=VCFv4.1\n')
                    chr2_file.write(date+'\n')
                    chr2_file.write('##contig=<ID='+str(j)+'>\n')
                    chr2_file.write('##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n')
                    chr2_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    chr2_file.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
                    chr2_file.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">\n')
                    chr2_file.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">\n')
                    chr2_file.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
                    header=('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
                    chr2_file.write(header+kk+'\n')
                    
                    chr_sites=np.where(chr_ord==j) #keeps snps in this chr for all inds
                    geno_pos=usable_pos[chr_sites] # keep genomic positions of sites in chr1 
                    phase_chr=phase_arr[(chr_sites)][:,i1:] #keep phased genotype info of sites in chr1
                    chr_genos=g1[chr_sites][:,i1:]
                    chr_genos_tot=g1[chr_sites]
                    cov_arr=np.zeros((len(phase_chr), len_pop2_hap), dtype='object') #must be the hap size of pop2
                    a1_dp=np.zeros((len(phase_chr), len_pop2_hap), dtype='object')
                    a2_dp=np.zeros((len(phase_chr), len_pop2_hap), dtype='object')
                    


                    for inx, g in enumerate(phase_chr[:,indx-i1]): # all the snps for this individual 
                       # print(inx,g, 'person:',indx-i1, indx, 'chr:',j)


                        cov=np.random.poisson(float(m)) #coverage for this site

                        if g=='0|0':
                            a1=cov     
                            a2=0
                            flip1_p=np.random.uniform(0, 1, (a1))
                            flip_a1=len(np.where(flip1_p>1-0.010)[0]) # 1% chance of flip error. How many reads will flip
                            flip_a2=0
                            a1_err=a1-flip_a1
                            a2_err=flip_a2+flip_a1
            
                            cov_arr[inx, indx-i1]=cov
                            a1_dp[inx, indx-i1]=a1_err
                            a2_dp[inx, indx-i1]=a2_err
            
            
                        elif g=='1|1':
                            a1=0
                            a2=cov
                            flip2_p=np.random.uniform(0, 1,(a2))
                            flip_a2=len(np.where(flip2_p>1-0.010)[0]) # 1 % chance of flip error
                            flip_a1=0
                            a1_err=flip_a1+flip_a2
                            a2_err=a2-flip_a2
            
                            cov_arr[inx, indx-i1]=cov
                            a1_dp[inx,indx-i1]=a1_err
                            a2_dp[inx,indx-i1]=a2_err
                
            
                        elif g=='1|0' or g=='0|1':
                            a1=np.random.binomial(cov, 0.5)
                            a2=cov-a1
                            flip1_p=np.random.uniform(0, 1, (a1))
                            flip2_p=np.random.uniform(0, 1, (a2))
                            flip_a1=len(np.where(flip1_p>1-0.010)[0]) # 1 % chance of flip error
                            flip_a2=len(np.where(flip2_p>1-0.010)[0]) # 1 % chance of flip error
                            a1_err=a1-flip_a1+flip_a2
                            a2_err=a2-flip_a2+flip_a1    

                            cov_arr[inx,indx-i1]=cov
                            a1_dp[inx,indx-i1]=a1_err
                            a2_dp[inx,indx-i1]=a2_err        
            
            
                        pos=str(geno_pos[inx])
                        ref=str(np.unique(chr_genos_tot[inx,:])[0])

                    
                        if len(np.unique(chr_genos_tot[inx,:]))==2:
                            alt=str(np.unique(chr_genos_tot[inx,:])[1]) #Keep the alt allele


                    # For GL estimation
                        if (phase_chr[inx,indx-i1])[0]=='0':
                            a1=ref
                        elif (phase_chr[inx,indx-i1])[0]=='1':
                            a1=alt
                    
                        if (phase_chr[inx,indx-i1])[2]=='0':
                            a2=ref

                        elif (phase_chr[inx,indx-i1])[2]=='1':
                            a2=alt




                        # Simulate Qscores based on Viking data VK350 

                        #Simulate position od this site in the reads
                        read_pos=int(np.random.uniform(1,83))

                        # According to this position, simulate Quality score
                        if read_pos <= 5:
                            qual=np.random.poisson(33)

                        elif 5 < read_pos <= 19:
                            qual=np.random.poisson(36)

                        elif 19 < read_pos <= 79:
                            qual=np.random.poisson(37)

                        elif 79 < read_pos < 82:
                            qual=np.random.poisson(36)

                        
                    
                        if a2_dp[inx,indx-i1]==0 and a1_dp[inx,indx-i1]!=0: # if all the cov goes to one allele
                            gg1=[a1,qual]
                            d1=a1_dp[inx,indx-i1]
                            d2=a2_dp[inx,indx-i1]
                            log_GL=geno_caller_3GT([gg1 for i in range(d1)], ref, alt, all_dic)

                            if a1==ref and a2==alt:
                            	dp_ref=d1
                            	dp_alt=0
                            elif a1==alt and a2==ref: #1/0 all coverage goes to the alt allele
                                dp_ref=0
                                dp_alt=d1
                            elif a1==ref and a2==ref:
                            	dp_ref=d1
                            	dp_alt=0
                            elif a1==alt and a2==alt:
                            	dp_ref=0
                            	dp_alt=d1

                        
                        elif a1_dp[inx,indx-i1]==0 and a2_dp[inx,indx-i1]!=0: #if all the cov goes to the other allele

                            gg2=[a2,qual]
                            d1=a2_dp[inx,indx-i1]
                            d2=a2_dp[inx,indx-i1]
                            log_GL=geno_caller_3GT([gg2 for i in range(d2)], ref, alt, all_dic)

                            if a1==ref and a2==alt: #0/1
                                dp_ref=0
                                dp_alt=d2
                            elif a1==alt and a2==ref: #1/0
                            	dp_ref=d2
                            	dp_alt=0
                            elif a1==ref and a2==ref: #0/0
                            	dp_ref=d2
                            	dp_alt=0
                            elif a1==alt and a2==alt: #1/1
                            	dp_ref=0
                            	dp_alt=d2

                        elif a1_dp[inx,indx-i1]!=0 and a2_dp[inx,indx-i1]!=0:
                        	gg1=[a1,qual]
                        	d1=a1_dp[inx,indx-i1]
                        	gg2=[a2,qual]
                        	d2=a2_dp[inx,indx-i1]
                        	lista1=[gg1 for i in range(d1)]
                        	lista2=[gg2 for i in range(d2)]

                        	if a1==alt and a2==ref:
                        		log_GL=geno_caller_3GT(lista1+lista2, a2, a1, all_dic)
                        		dp_ref=d2
                        		dp_alt=d1
                        	elif a1==ref and a2==alt:
                        		log_GL=geno_caller_3GT(lista1+lista2, a1, a2, all_dic)
                        		dp_ref=d1
                        		dp_alt=d2
                        	elif a1==ref and a2==ref:
                        		a2=alt
                        		gg1=[a1,qual]
                        		gg2=[a2,qual]
                        		d1=a1_dp[inx,indx-i1]
                        		d2=a2_dp[inx,indx-i1]
                        		dp_ref=d1
                        		dp_alt=d2
                        		lista1=[gg1 for i in range(d1)]
                        		lista2=[gg2 for i in range(d2)]
                        		log_GL=geno_caller_3GT(lista1+lista2,a1,a2,all_dic)
                        	elif a1==alt and a2==alt: #1/1 both with coverage (due to of seq error)
                        		a1=ref 
                        		gg1=[a1,qual]
                        		gg2=[a2,qual]
                        		d1=a1_dp[inx,indx-i1]
                        		d2=a2_dp[inx,indx-i1]
                        		dp_ref=d1
                        		dp_alt=d2
                        		lista1=[gg1 for i in range(d1)]
                        		lista2=[gg2 for i in range(d2)]
                        		log_GL=geno_caller_3GT(lista1+lista2,a1,a2,all_dic)

                       
                                                
                        if cov_arr[inx,indx-i1]!=0:

                            raw_PL=[math.log10(10**i)*-10 for i in log_GL] # Raw Genotype likelihood
                            norm_PL=[(round(i-min(raw_PL))) for i in raw_PL] # Normalized Genotype Likelihoods 
                            GQ=int(sorted(norm_PL)[1])
                            norm_PL=str(norm_PL).replace('.0', '')
                            norm_PL=norm_PL.replace(' ', '')

                        elif cov_arr[inx,indx-i1]==0: # If coverage is 0 (missing site)
                            phase_chr[inx,indx-i1]='./.'
                            GQ='.'
                            cov_arr[inx,indx-i1]='0'
                            a1_dp[inx,indx-i1]='0'
                            a2_dp[inx,indx-i1]='0'
                            dp_ref='0'
                            dp_alt='0'
                            norm_PL='0,0,0' 

                        # WRITE VCF FILE LINES   

                        if len(np.unique(chr_genos_tot[inx,:]))==2:
                            alt=str(np.unique(chr_genos_tot[inx,:])[1])
                            out= str(j)+ '\t'+pos+'\t'+'snp_'+str(j)+pos+'\t'+ref+'\t'+alt+'\t.\t.\tPR\t'+'GT:GQ:DP:AD:PL'+'\t'+str(str(phase_chr[inx,indx-i1]).split())+':'+str(GQ)+':'+str(cov_arr[inx,indx-i1])+':'+str(dp_ref)+','+str(dp_alt)+':'+norm_PL+'\n'
                            out=out.replace('[','')
                            out=out.replace(']','')
                            out=out.replace('"','')
                            out=out.replace("'",'')
                            out=out.replace("|",'/')
                            out=out.replace(" ",'\t')
                            out=out.replace(".0", '')
                            chr2_file.write(out)
                       

print datetime.datetime.now() - startTime 

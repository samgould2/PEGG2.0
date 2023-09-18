#!/usr/bin/env python
# coding: utf-8

# In[8]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
import gzip
from Bio import AlignIO
import Bio.Align
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import matplotlib.patches as patches
import re
import seaborn as sns
#import matplotlib.patches as patches
from matplotlib.patches import Polygon
from cyvcf2 import VCF
pd.set_option('display.max_columns', 50)


#-----------functions------------
def genome_loader(filepath_gz):
    """
    Takes in filepath of human genome (GrCH37 or GrCh38) and returns records and index_list for PEGG parsing.
    
    Parameters
    -----------
    filepath_gz
        *type = str*
        
        The filepath to the .gz file holding the reference genome file.
    
    """
    #------loading in reference genome and organizing it into a 2-d list by chromosome---------------------

    with gzip.open(filepath_gz, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta")) #about 4 Gb in  memory
        #records = list that contains sequences split up by chromosome (and intrachromosome splits up to some size)

    #filtering out alternative sequences to only select consensus chromosome sequences
    wrong = ["alternate", "unplaced", "unlocalized", "patch"]
    badlist = []
    for key in wrong:
        for i in records:
            ii = i.description
            if key in ii:
                badlist.append(ii)

    #creating an 
    filtered = []
    index_list = []
    for idx, i in enumerate(records):
        ii = i.description
        if ii not in badlist:
            filtered.append(ii)
            index_list.append(idx)
            
    return records, index_list


def PAM_finder(mutant_input, PAM, RTT_length, mut_idx, records, index_list):
    """Identifies the location of PAM sequences on the + and - strand.
    Returns a 2-d array containing marked locations of PAM sequence start locations on + and - strand.
    
    Parameters
    ----------
    mutant_input
        *type = pd.DataFrame*
        
        A dataframe containing the input mutations from which selections are made to generate pegRNAs.
        See documentation for precise qualities of this dataframe
        
    PAM
        *type = str*
        
        String designating what PAM sequence to search for. Not all PAM sequences exhaustively tested.
        Normally this will be equal to 'NGG' for most existing prime editors.
        
    RTT_length
        *type = int*
        
        Specifies the reverse transcriptase template (RTT) length of the pegRNA being generated.
        
    mut_idx
        *type = int*
        
        Index of mutation within mutant_input. Defines which mutation to identify PAM sequences for.
        
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records (reference genome). See documentation for precise requirements.
    
    """
    #----------------Loading in mutation information & sequence -----------------#
    mut = mutant_input.iloc[[mut_idx]]

    mut_type = mut['Variant_Type'] 

    #getting gene and start/end position of mutation
    gene = mut['Hugo_Symbol'].values[0]

    s = mut['Start_Position'].values[0]
    e = mut['End_Position'].values[0]
    size_mut = (e-s)+1 #size of mutation; not correct for insertions though

    ##loading in genome information
    ref_allele = mut['Reference_Allele'].values[0]
    mut_allele = mut['Tumor_Seq_Allele2'].values[0]

    seq_start = mut['Start_Position'].values[0]
    seq_end = mut['End_Position'].values[0]
    chromosome = mut['Chromosome'].values[0]
    
    if chromosome=='X':
        chrom = 22
    else:
        chrom = int(chromosome)-1
        
    seq1 = records[index_list[chrom]].seq
    #seq1 = records[index_list[int(chromosome)-1]].seq

    #---------------Loading in sequences for PAM Searching------------------#
    
    search_size = RTT_length - len(PAM)-1 #need to modify this for insertion/deletions...
    #size mut doesn't capture the size of insertions, only the size of deletions

    plus_search = seq1[seq_start-1-search_size : seq_end+search_size].upper()
    minus_search = plus_search.complement().upper()
    
    mut_start_idx = 1+search_size
    mut_end_idx = 1+search_size+size_mut #not accurate for insertions; does work for indexing though...

    plus_search1 = plus_search[:mut_start_idx+3+len(PAM)-1]
    minus_search1 = minus_search[mut_start_idx-3-len(PAM):]

    #---------------PAM Searching------------------#

    #replacing N with regex symbol
    PAM_regex = PAM.replace('N', '/*.')

    PAM_search_plus = re.compile('(?=(' + PAM_regex + '))', re.IGNORECASE)

    iterator_plus = PAM_search_plus.finditer(str(plus_search1))
    PAM_starts_plus = [match.start() for match in iterator_plus]


    PAM_minus = PAM[::-1]#reversing it
    PAM_regex_minus = PAM_minus.replace('N', '/*.')

    PAM_search_minus = re.compile('(?=(' + PAM_regex_minus + '))', re.IGNORECASE)

    iterator_minus = PAM_search_minus.finditer(str(minus_search1))
    PAM_starts_minus = [match.start() for match in iterator_minus]
    #since things are flipped on minus strand, adding len(PAM) to get the true "start" to the PAM

    PAM_starts_minus = np.asarray(PAM_starts_minus) + len(PAM) + (mut_start_idx-3-len(PAM))#and correct for indexing

    return np.asarray([np.asarray(PAM_starts_plus), PAM_starts_minus], dtype='object')-mut_start_idx

    #return index of START of PAM sequence relative to the start of the mutation

def PAM_finder_multi_index(mutant_input,PAM, RTT_length, mut_idx_list, records, index_list):
    """Identifies the location of PAM sequences on the + and - strand for a LIST OF INPUT MUTATIONS. 
    Returns a list of 2-d arrays containing marked locations of PAM sequence start locations on + and - strand.
    
    Parameters
    ----------
    mutant_input
        *type = pd.DataFrame*
        
        A dataframe containing the input mutations from which selections are made to generate pegRNAs.
        See documentation for precise qualities of this dataframe
        
    PAM
        *type = str*
        
        String designating what PAM sequence to search for. Not all PAM sequences exhaustively tested.
        Normally this will be equal to 'NGG' for most existing prime editors.
        
    RTT_length
        *type = int*
        
        Specifies the reverse transcriptase template (RTT) length of the pegRNA being generated.
        
    mut_idx_list
        *type = list*
        
        List of indeces of mutations within mutant_input that user wants to identify PAM sequences for.
        
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records (reference genome). See documentation for precise requirements.
    
    """
    list_of_marked_arrays = []
    
    for mut_idx in mut_idx_list:
        marked_arrays = PAM_finder(mutant_input,PAM, RTT_length, mut_idx, records, index_list)
        
        list_of_marked_arrays.append(marked_arrays)
    
    return list_of_marked_arrays


def target_design_w_flank(PAM_location, seq_start, seq_end, seq1, max_size, strand):
    """Generates a synthetic version of the endogenous target site ("sensor" site). 
    This sensor region can be used as a proxy readout of editing outcomes at the endogenous locus if 
    included in the oligo with the pegRNA. For more information about sensor sites see: 
    https://www.nature.com/articles/s41587-021-01172-3
    
    Parameters
    ----------
    PAM_location
        *type = int*
        
        Int of location of PAM sequence (generated from PAM_finder function).
        
    seq_start
        *type = int*
        
        Start site of mutation of interest.
        
    seq_end
        *type = int*
        
        End site of mutation of interest.
        
    seq1
        *type = str*
        
        Chromosome sequence that corresponds to location of mutation of interest  .   
        
    max_size
        *type = int*
        
        Size of synthetic target site. Reccomended size = 60 nt. 
        
    strand
        *type = str*
        
        Strand of PAM sequence. Options = '+' or '-'.
    """
    
    #automated target design, but also includes flanking sequence (if allowed by max size...)
    if strand == '+':
        protospacer_w_PAM = seq1[PAM_location-20:PAM_location+3].upper()
        
        #calc distance between last G of NGG in PAM sequence and mutation
        PAM_mut_dist = seq_end - (PAM_location+2)
        
        if PAM_mut_dist>0: #if it's ahead of the PAM sequence
            incl_mut_seq = seq1[PAM_location+3:PAM_location+2+PAM_mut_dist].upper()
            
            target_sequence = str(protospacer_w_PAM) +  str(incl_mut_seq)
            
        else:
            target_sequence = str(protospacer_w_PAM)
    
    else:
        #input = complement for this
        protospacer_w_PAM = seq1[PAM_location-20:PAM_location+3].upper()
        
        #calc distance between last G of NGG in PAM sequence and mutation
        PAM_mut_dist = seq_start - (PAM_location+2) + 1 #start instead of end because of flipped orientation
        
        if PAM_mut_dist>0: #if it's ahead of the PAM sequence
            incl_mut_seq = seq1[PAM_location+3:PAM_location+2+PAM_mut_dist].upper()
            
            target_sequence = str(protospacer_w_PAM) +  str(incl_mut_seq)
            
        else:
            target_sequence = str(protospacer_w_PAM)

        s = Bio.Seq.Seq(target_sequence).reverse_complement()
    
        target_sequence = str(s)
        
    target_len = len(target_sequence)
    
    #––––––––––––now adding the flanking sequences–––––––––––––––––––––––––––––––––––#
    left_over = max_size-target_len
    
    if left_over>0:
        
        if left_over>5:
            #add 5 to downstream of PAM and the rest upstream
            if strand == '+':
                
                pam_side_l = seq1[PAM_location-25:PAM_location-20].upper()
                if PAM_mut_dist>0:
                    mut_side_r = seq1[PAM_location+2+PAM_mut_dist:PAM_location+2+PAM_mut_dist+(left_over-5)]
                elif PAM_mut_dist<=0:
                    mut_side_r = seq1[PAM_location+3:PAM_location+3+(left_over-5)]
                                      
                target_sequence = str(pam_side_l) + str(target_sequence) + str(mut_side_r)
                
                
                
                
            if strand == '-': #minus strand
                
                pam_side_l = seq1[PAM_location-25:PAM_location-20].upper()
                if PAM_mut_dist>0:
                    mut_side_r = seq1[PAM_location+2+PAM_mut_dist:PAM_location+2+PAM_mut_dist+(left_over-5)]
                elif PAM_mut_dist<=0:
                    mut_side_r = seq1[PAM_location+3:PAM_location+3+(left_over-5)]
                
                target_sequence = str(mut_side_r.reverse_complement()) + str(target_sequence) + str(pam_side_l.reverse_complement())


                
                #s = Bio.Seq.Seq(target_sequence).reverse_complement()
    
                #target_sequence = str(s)
                
                
                
            
        else:
            #only add to upstream of PAM (near edit)
            if strand == '+':
                
                if PAM_mut_dist>0:
                    mut_side_r = seq1[PAM_location+2+PAM_mut_dist:PAM_location+2+PAM_mut_dist+(left_over)]
                else:
                    mut_side_r = seq1[PAM_location+3:PAM_location+3+(left_over)]
                                      
                target_sequence = str(target_sequence) + str(mut_side_r)
                
                
                
            else: #minus strand
                
                if PAM_mut_dist>0:
                    mut_side_r = seq1[PAM_location+2+PAM_mut_dist:PAM_location+2+PAM_mut_dist+(left_over-5)]
                elif PAM_mut_dist<=0:
                    mut_side_r = seq1[PAM_location+3:PAM_location+3+(left_over-5)]
                
                target_sequence = str(mut_side_r.reverse_complement()) + str(target_sequence)


            
    else:
        target_sequence=target_sequence
            
            
        
        
    return target_sequence.upper()

def minus_seq_generator(records, index_list):
    
    """
     Simple function that generates reverse complement of chromosome sequences.
     Reduces run time significantly by preventing undertaking this action multiple times.
     Returns list corresponding to reverse complement of records.
     
    Parameters
    ----------
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records (reference genome). See documentation for precise requirements.
    
    """
    
    #taking reverse complement of chromosomes to get minus strand in 5' to 3' orientation
    minus_seqs = [] 
    for i in range(len(index_list)):
        seq1 = records[index_list[i]].seq
        seq2 = seq1.reverse_complement()
        minus_seqs.append(seq2)
    
    return minus_seqs

def mutation_consequences(mutant_input, mut_idx_list, records, index_list):
    """
    Calculates the consequences of mutations at their genomic locus.
    Returns a list of genomic regions 50nt up/downstream with desired edit included (~100 nt total size +/- size of INS or DEL).
    Can be used for quantifying editing outcomes and other alignment purposes.
    
    Parameters
    ----------
    mutant_input
        *type = pd.DataFrame*
        
        A dataframe containing the input mutations from which selections are made to generate pegRNAs.
        See documentation for precise qualities of this dataframe
        
    mut_idx_list
        *type = list*
        
        List of indeces of mutation within mutant_input that user wants to perform function on (i.e. a subset of the mutations).
        
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records (reference genome). See documentation for precise requirements.
    
    """
    
    
    output = [] #consequences of mutation at the genomic locus
    for idx, val in enumerate(mut_idx_list):

        mut = mutant_input.iloc[[val]]

        mut_type = mut['Variant_Type'].values[0]

        #getting gene and start/end position of mutation
        gene = mut['Hugo_Symbol'].values[0]

        s = mut['Start_Position'].values[0]
        e = mut['End_Position'].values[0]
        size_mut = (e-s)+1 #size of mutation; not correct for insertions though

            ##loading in genome information
        ref_allele = mut['Reference_Allele'].values[0]
        mut_allele = mut['Tumor_Seq_Allele2'].values[0]

        seq_start = mut['Start_Position'].values[0]
        seq_end = mut['End_Position'].values[0]
        chromosome = mut['Chromosome'].values[0]

        if chromosome=='X':
            chrom = 22
        else:
            chrom = int(chromosome)-1

        seq1 = records[index_list[chrom]].seq  
            
        if mut_type=='DEL':
            if mut_allele == '-':
                out = str(seq1[seq_start-50:seq_start-1]) +str(seq1[seq_end:seq_end+50])
            else:
                out = str(seq1[seq_start-50:seq_start-1]) + str(mut_allele) + str(seq1[seq_end:seq_end+50])
        
        elif mut_type=='INS':
            if ref_allele == '-':
                out = str(seq1[seq_start-50:seq_start-1]) + str(mut_allele) + str(seq1[seq_start-1:seq_start+50])
            else:
                out = str(seq1[seq_start-50:seq_start-1]) + str(mut_allele) + str(seq1[seq_end:seq_end+50])

            
        else: #SNPs and ONPs
            out = str(seq1[seq_start-50:seq_start-1]) + str(mut_allele) + str(seq1[seq_end:seq_end+50])
         
        output.append(out.upper())
    
        
    return output

def pegRNA_generator(mutant_input, PBS_length, RTT_length, PAM, marked_array_list, mut_idx, records, index_list, minus_seqs):    
    """Generates all possible pegRNAs matching the input specifications.
    Returns a dataframe containing the pegRNAs and some of their properties quantified.
    
    Parameters
    ----------
    mutant_input
        *type = pd.DataFrame*
        
        A dataframe containing the input mutations from which selections are made to generate pegRNAs.
        See documentation for precise qualities of this dataframe
     
    PBS_length
        *type = int*
        
        Specifies the primer binding sequence (PBS) length of the pegRNA being generated.
        
    RTT_length
        *type = int*
        
        Specifies the reverse transcriptase template (RTT) length of the pegRNA being generated. 
        Must match input into PAM_finder (so that marked_array_list matches).
    
    PAM
        *type = str*
        
        String designating what PAM sequence to search for. Not all PAM sequences exhaustively tested.
        Normally this will be equal to 'NGG' for most existing prime editors.
        
    marked_array_list
        *type = list*
        
        List containing locations of PAM sequences GENERATED BY PAM_finder function (and/or PAM_finder_multi_index)
    
    mut_idx
        *type = list or int*
        
        List of indeces (or single index (type=int)) of mutations within mutant_input that user wants to perform function on (i.e. a subset of the mutations).        
        
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records list (reference genome). See documentation for precise requirements.
        
    minus_seqs
        *type = list*
        
        List containing reverse complement of records (reference genome). Generated by minus_seq_generator function.
    
    """
    
    #start by initializing an ugly bunch of lists for storing information
    PAM_locs = []
    PAM_seqs = []
    PAM_strand = []
    protos = []
    PBS_seqs = []
    RT_seqs = []
    PBS_RT_5to3 = []
    PBS_len = []
    RT_len = []
    mut_to_3_distance = []
    nick_dist = []
    
    target_seq = []
    target_size = []

    impact_idx = []
    gene_name = []
    chrom_num = []
    start1 = []
    end1 = []
    ref1 = []
    mut1 = []
    mut_type1 = []
    
    #and creating the function to actually create pegRNAs...
    if type(mut_idx)==int:
        mut_idx  = [mut_idx]
    
    
    max_size = 60 #should make this an input
    
    for idx, val in enumerate(mut_idx):
        
        mut = mutant_input.iloc[[val]]
    
        mut_type = mut['Variant_Type'].values[0]

        #getting gene and start/end position of mutation
        gene = mut['Hugo_Symbol'].values[0]

        s = mut['Start_Position'].values[0]
        e = mut['End_Position'].values[0]
        size_mut = (e-s)+1 #size of mutation; not correct for insertions though

            ##loading in genome information
        ref_allele = mut['Reference_Allele'].values[0]
        mut_allele = mut['Tumor_Seq_Allele2'].values[0]

        seq_start = mut['Start_Position'].values[0]
        seq_end = mut['End_Position'].values[0]
        chromosome = mut['Chromosome'].values[0]
        
        if chromosome=='X':
            chrom = 22
        else:
            chrom = int(chromosome)-1
        
        seq1 = records[index_list[chrom]].seq

        ##----------Generate pegRNAs--------------##
        
        
        

        #start with plus strand
        for i in marked_array_list[idx][0]:
            chrom_num.append('chr' + str(chromosome))
            
            PAM_location = seq_start+i
            
            strand = '+'
            target_sequence = target_design_w_flank(PAM_location, seq_start, seq_end, seq1, max_size, strand)
            target_seq.append(target_sequence)
            target_size.append(len(target_sequence))
            
            PAM_sequence = seq1[PAM_location:PAM_location+len(PAM)].upper()

            protospacer = seq1[PAM_location-20:PAM_location].upper() #5' to 3'
            PBS = seq1[PAM_location-4-PBS_length+1:PAM_location-3].complement().upper() #3' to 5'


            ###NEED TO TAKE INTO ACCOUNT VARIANT TYPE FOR RT TEMPLATE DESIGN
            if mut_type=='SNP': #for now just worry about SNPs
                RT = str(seq1[PAM_location-3:seq_start-1])+str(mut_allele)+str(seq1[seq_start:PAM_location-3+RTT_length])
                RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                mut_to_3 = abs((PAM_location-3+(RTT_length))-seq_start)

            elif mut_type=='DEL':

                ss = abs((PAM_location-3) - (seq_start-1))

                if mut_allele =='-':
                    del_size = abs(seq_end-seq_start)+1

                    RT_no_end = str(seq1[PAM_location-3:seq_start-1])

                    if len(RT_no_end)>RTT_length:
                        RT = 'RTT not long enough'
                    else:
                        RT = str(seq1[PAM_location-3:seq_start])+str(seq1[seq_start+del_size:PAM_location-3+(RTT_length+del_size)])
                        RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                        mut_to_3 = abs((PAM_location-3+(RTT_length+del_size))-seq_start)-del_size



                else:

                    del_size = abs(seq_end-seq_start)-len(mut_allele)+1

                    RT_no_end = str(seq1[PAM_location-3:seq_start-1])+str(mut_allele)

                    if len(RT_no_end)>RTT_length:
                        RT = 'RTT not long enough'
                    else:
                        RT = str(seq1[PAM_location-3:seq_start-1])+str(mut_allele)+str(seq1[seq_start+del_size:PAM_location-3+(RTT_length+del_size)])
                        RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                        mut_to_3 = abs((PAM_location-3+(RTT_length+del_size))-seq_start)-del_size

            elif mut_type =='INS':

                if ref_allele =='-':
                    ins_size = len(mut_allele)

                    RT_no_end = str(seq1[PAM_location-3:seq_start-1])+str(mut_allele)

                    if len(RT_no_end)>RTT_length:
                        RT = 'RTT not long enough'
                    else:
                        RT = str(seq1[PAM_location-3:seq_start-1])+str(mut_allele)+str(seq1[seq_start-1:PAM_location-3+(RTT_length-ins_size)])
                        RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                else:

                    ins_size = len(mut_allele)-len(ref_allele)

                    RT_no_end = str(seq1[PAM_location-3:seq_start-1])+str(mut_allele)

                    if len(RT_no_end)>RTT_length:
                        RT = 'RTT not long enough'
                    else:
                        RT = str(seq1[PAM_location-3:seq_start-1])+str(mut_allele)+str(seq1[seq_start+len(ref_allele)-1:PAM_location-3+(RTT_length-ins_size)])
                        RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                mut_to_3 = abs((PAM_location-3+(RTT_length-ins_size))-seq_start)

            elif mut_type ==('ONP' or 'DNP'):
                
                
                del_size = abs(seq_end-seq_start)-len(mut_allele)+1

                RT_no_end = str(seq1[PAM_location-3:seq_start-1])+str(mut_allele)

                if len(RT_no_end)>RTT_length:
                    RT = 'RTT not long enough'
                else:
                    RT = str(seq1[PAM_location-3:seq_start-1])+str(mut_allele)+str(seq1[seq_end+del_size:PAM_location-3+(RTT_length+del_size)])
                    RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                mut_to_3 = abs((PAM_location-3+(RTT_length+del_size))-seq_start)-del_size

            #elif mut_type == 'DNP':

            else:
                continue


            #combining PBS and RT into a single DNA sequence in the correct orientation
            if RT== 'RTT not long enough':
                PBS_RT = 'n/a'
            else: 
                PBS_RT = (PBS+RT).reverse_complement().complement() #given in 5' to 3' orientiation

                
            #distance to nick
            nd = RTT_length-mut_to_3
                
            #pegRNA info
            PAM_locs.append(PAM_location)
            PAM_strand.append('+')
            PAM_seqs.append(str(PAM_sequence))
            protos.append(str(protospacer)) #5' to 3'
            PBS_seqs.append(str(PBS)) #3' to 5'
            PBS_len.append(PBS_length)
            RT_seqs.append(str(RT)) #3' to 5'
            RT_len.append(len(RT))
            PBS_RT_5to3.append(str(PBS_RT))
            mut_to_3_distance.append(mut_to_3)
            nick_dist.append(nd)
            
            #mutation info
            ref1.append(ref_allele)
            mut1.append(mut_allele)
            mut_type1.append(mut_type)
            gene_name.append(gene)
            impact_idx.append(val)
            start1.append(seq_start)
            end1.append(seq_end)

        #then do the min    
        for i in marked_array_list[idx][1]:
            chrom_num.append('chr' + str(chromosome))
            
            #taking ther reverse complement of everything and then flipping it back to the correct orientation
            seq_start = mut['Start_Position'].values[0]
            seq_end = mut['End_Position'].values[0]

            mut_allele = mut['Tumor_Seq_Allele2'].values[0]

            mut_allele2 = Bio.Seq.Seq(mut_allele).reverse_complement()
            
            seq_start = len(seq1)-seq_start
            seq_end = len(seq1)-seq_end
            
            PAM_location = seq_start-i #minus instead of plus because of flipped orientation
            
            seq2 = minus_seqs[chrom]
            
            #looking for target
            strand = '-'
            target_sequence = target_design_w_flank(PAM_location, seq_start, seq_end, seq2, max_size, strand)
            target_seq.append(target_sequence)
            target_size.append(len(target_sequence))
            
            PAM_sequence = seq2[PAM_location:PAM_location+len(PAM)].upper()

            protospacer = seq2[PAM_location-20:PAM_location].upper() #5' to 3'
            
            PBS = seq2[PAM_location-4-PBS_length+1:PAM_location-3].complement().upper() #3' to 5'

            
            ##
            #need to flip indexing to get things to work__________________________
            seq_start1 = seq_end
            seq_end1 = seq_start
                
            seq_start = seq_start1
            seq_end=seq_end1

            ###NEED TO TAKE INTO ACCOUNT VARIANT TYPE FOR RT TEMPLATE DESIGN
            if mut_type=='SNP': #for now just worry about SNPs
                RT = str(seq2[PAM_location-3:seq_start])+str(mut_allele2)+str(seq2[seq_start+1:PAM_location-3+RTT_length])
                RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                mut_to_3 = abs((PAM_location-3+(RTT_length))-seq_start)

            elif mut_type=='DEL':
                

                ss = abs((PAM_location-3) - (seq_start-1))

                if mut_allele =='-':
                    del_size = abs(seq_end-seq_start)+1

                    RT_no_end = str(seq2[PAM_location-3:seq_start-1])

                    if len(RT_no_end)>RTT_length:
                        RT = 'RTT not long enough'
                    else:
                        RT = str(seq2[PAM_location-3:seq_start])+str(seq2[seq_start+del_size:PAM_location-3+(RTT_length+del_size)])
                        RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                        mut_to_3 = abs((PAM_location-3+(RTT_length+del_size))-seq_start)-del_size



                else:
                    
                    del_size = abs(seq_end-seq_start)-len(mut_allele)+1

                    RT_no_end = str(seq2[PAM_location-3:seq_start])+str(mut_allele2)

                    if len(RT_no_end)>RTT_length:
                        RT = 'RTT not long enough'
                        
                    else:
                        #print(Bio.Seq.Seq(seq2[PAM_location-3:seq_start]).reverse_complement())
                        #print(mut_allele2.reverse_complement())
                        RT = str(seq2[PAM_location-3:seq_start])+str(mut_allele2)+str(seq2[seq_start+1+del_size:PAM_location-3+(RTT_length+del_size)])
            
                        RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                        mut_to_3 = abs((PAM_location-3+(RTT_length+del_size))-seq_start)-del_size

            elif mut_type =='INS':

                if ref_allele =='-':
                    ins_size = len(mut_allele)

                    RT_no_end = str(seq2[PAM_location-3:seq_start-1])+str(mut_allele2)

                    if len(RT_no_end)>RTT_length:
                        RT = 'RTT not long enough'
                    else:
                        RT = str(seq2[PAM_location-3:seq_end+1])+str(mut_allele2)+str(seq2[seq_end+1:PAM_location-3+(RTT_length-ins_size)])
                        RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                else:

                    ins_size = len(mut_allele)-len(ref_allele)

                    RT_no_end = str(seq2[PAM_location-3:seq_start-1])+str(mut_allele2)

                    if len(RT_no_end)>RTT_length:
                        RT = 'RTT not long enough'
                    else:
                        RT = str(seq2[PAM_location-3:seq_start])+str(mut_allele2)+str(seq2[seq_start+len(ref_allele):PAM_location-3+(RTT_length-ins_size)])
                        RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                mut_to_3 = abs((PAM_location-3+(RTT_length-ins_size))-seq_start)

            elif mut_type ==('ONP' or 'DNP'): #trying to get this to work...
                
                
                del_size = abs(seq_end-seq_start)-len(mut_allele)+1

                RT_no_end = str(seq2[PAM_location-3:seq_start])+str(mut_allele2)

                if len(RT_no_end)>RTT_length:
                    RT = 'RTT not long enough'
                    
                else:
                    #print(Bio.Seq.Seq(seq2[PAM_location-3:seq_start]).reverse_complement())
                    #print(mut_allele2.reverse_complement())

                    RT = str(seq2[PAM_location-3:seq_start])+str(mut_allele2)+str(seq2[seq_end+1+del_size:PAM_location-3+(RTT_length+del_size)])

                    RT = Bio.Seq.Seq(RT).complement().upper() #3' to 5'

                mut_to_3 = abs((PAM_location-3+(RTT_length+del_size))-seq_start)-del_size

            #elif mut_type == 'DNP':

            else:
                continue
            

            #combining PBS and RT into a single DNA sequence in the correct orientation
            if RT== 'RTT not long enough':
                PBS_RT = 'n/a'
            else: 
                RT = RT.reverse_complement().complement()
                PBS = PBS.reverse_complement().complement()
                
                PBS_RT = (RT+PBS) #given in 5' to 3' orientiation


            
            
            nd = RTT_length-mut_to_3
            
            seq_start = mut['Start_Position'].values[0]
            seq_end = mut['End_Position'].values[0]
            PAM_location = seq_start+i
            
            #pegRNA info
            PAM_locs.append(PAM_location)
            PAM_strand.append('-')
            PAM_seqs.append(str(PAM_sequence))
            protos.append(str(protospacer)) #5' to 3'
            PBS_seqs.append(str(PBS)) #3' to 5'
            PBS_len.append(PBS_length)
            RT_seqs.append(str(RT)) #3' to 5'
            RT_len.append(len(RT))
            PBS_RT_5to3.append(str(PBS_RT))
            mut_to_3_distance.append(mut_to_3)
            nick_dist.append(nd)

            
            #mutation info
            ref1.append(ref_allele)
            mut1.append(mut_allele)
            mut_type1.append(mut_type)
            gene_name.append(gene)
            impact_idx.append(val)
            start1.append(seq_start)
            end1.append(seq_end)
            

            
            
    list_of_tuples = list(zip(impact_idx,gene_name, chrom_num, start1, end1, mut_type1, ref1, mut1,
                             PAM_locs, PAM_strand, PAM_seqs, protos, PBS_seqs,PBS_len, RT_seqs, RT_len, PBS_RT_5to3,
                             mut_to_3_distance,nick_dist, target_seq, target_size))

    pegRNA_df = pd.DataFrame(list_of_tuples, columns = ['mutant index', 'gene', 'chr','start', 'end', 'variant_type', 'ref_allele',
                                           'mut_allele', 'PAM location', 'PAM strand','PAM', 'protospacer', 'PBS',
                                            'PBS length', 'RTT', 'RTT length','PBS_RTT_5to3', "distance mut to 5' RTT","distance to nick", 
                                                        "target sequence", "target length"] )
    
    #filter out pegRNAs with too long of an RTT
    pegRNA_df = pegRNA_df[pegRNA_df["RTT"] != 'RTT not long enough']
    pegRNA_df = pegRNA_df[pegRNA_df["RTT length"] == RTT_length]


    return pegRNA_df.reset_index(drop=True)


#------------ranking and filtration functions--------
def on_off_target(pegRNA_df, chrom_dict):
    """Calculates the on/off-target scores for protospacers occuring in exons.
    If chrom_dict is not provided, set it equal to 'none' and ranking/filtration will proceed without on/off target score.
    
    Parameters
    ----------
    pegRNA_df
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()
     
    chrom_dict
        *type = dict*
        
        A dictionary containing the on/off-target scores for each of the exonic regions in each chromosome.
        Keys = chromosome identifiers. See documentation for more information and access to provided reference file.
        
        If no reference file available SET THIS EQUAL TO 'none' and ranking will proceed without on/off target score.
    
    
    """
    
    if chrom_dict != 'none':
        MIT_score = []
        CFD = []

        for idx, i in pegRNA_df.iterrows():
            chrom = i['chr']
            pam_site = i['PAM location']
            strand = i['PAM strand']

            df = chrom_dict[chrom]

            if strand=='+':
                row = df[df['thickEnd']==pam_site]

            elif strand=='-':
                row = df[df['thickStart']==pam_site]


            mit = row['score']
            fusi = row['fusi']
            if len(mit)==0:
                mit='none'
                fusi='none'

                MIT_score.append(mit)
                CFD.append(fusi)

            else:
                #fusi =fusi.values[0]
                mit=mit.values[0]
                fusi = fusi.values[0]

                MIT_score.append(mit)
                CFD.append(int(fusi.split('%')[0]))
     
    if chrom_dict=='none':
        MIT_score = ['none']*len(pegRNA_df)
        CFD = ['none']*len(pegRNA_df)
    
    pegRNA_df['MIT specificity']=MIT_score
    pegRNA_df['CFD efficiency percentile']=CFD #doench 2016 score
    
    return pegRNA_df
        
    
def other_filtration(pegRNA_df):

    """
    Classifies properties of pegRNAs in pegRNA_df besides on/off target score. Returns pegRNA_df with properties added.
    These properties are checking for:
    
    (1) u6 terminator (polyT sequence)
    
    (2) PBS GC content
    
    (3) RE site presence -- EcoRI and Esp3I sites (used for cloning in the provided adapter sequences)
    
    (4) Checks whether the last templated nucleotide is C (reccomended against in original Anzalone et al. paper)
    
    Parameters
    ----------
    pegRNA_df
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()
    """
    contains_terminator = []
    PBS_GC_content = []
    RE_site = []
    last_templated_c = []
    
    #making list of prohibited RE sites
    Esp3I = 'CGTCTC'
    Esp3I_comp = str(Bio.Seq.Seq(Esp3I).reverse_complement())
    EcoRI = 'GAATTC'
    EcoRI_comp = str(Bio.Seq.Seq(EcoRI).reverse_complement())
    RE_sites = [Esp3I, Esp3I_comp, EcoRI, EcoRI_comp]
    
    #check for u6 terminator sequence "TTTT" in pbs or protospacer
    for idx, i in pegRNA_df.iterrows():
        pbs_rtt = i["PBS_RTT_5to3"]
        proto = i["protospacer"]
        pbs = i["PBS"]
        target = i['target sequence']
        
        terminator_sequence = 'TTTT'
        
        if terminator_sequence in (pbs_rtt or proto):
            contains_terminator.append('yes')
        else:
            contains_terminator.append('no')
            
            
        #and calculate PBS GC content
        total_gc = pbs.count('G')+pbs.count('C')
        perc = total_gc/len(pbs)
        PBS_GC_content.append(perc)
        
        #finally check for restriction enzyme sites

        res_proto = any(re_site in proto for re_site in RE_sites)
        res_pbs_rtt = any(re_site in pbs_rtt for re_site in RE_sites)
        res_target = any(re_site in target for re_site in RE_sites)

        if (res_proto or res_pbs_rtt or res_target)==True:
            RE_site.append('yes')
        else:
            RE_site.append('no')    
        
        
        #lastly look at first 3' extension base
        if pbs_rtt[0].upper()=='C':
            last_templated_c.append('yes')
        else:
            last_templated_c.append('no')
            
    
    pegRNA_df['contains terminator?']=contains_terminator
    pegRNA_df['contains RE site?']=RE_site
    pegRNA_df['last templated base C?']= last_templated_c
    pegRNA_df['PBS GC content']=PBS_GC_content
    
    
    return pegRNA_df


def composite_score(pegRNA_df, RE_filter=True):
    """
    Calculates composite score of pegRNA properties for ranking purposes.
    This composite score is a weighted linear sum. Weights will be updated when training data is provided in future versions of PEGG.
    Requires that other_filtration() and on_off_target() are run on pegRNA_df first
    Returns pegRNA_df with composite score added added.

    Parameters
    ----------
    pegRNA_df
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()

    RE_filter
        *type = Boolean (True or False)

        If = True, filter for Esp3I and EcoRI sites. If = False, do not include this in the composite score...
    """
    
    
    x2 = list(pegRNA_df['MIT specificity'])
    if 'none' in x2: #checking for guides that aren't in region
        
        
        #numerical factors
        x1 = np.array(pegRNA_df["distance mut to 5' RTT"]) #want this to be large (edit close to nick site)
        x2 = list(pegRNA_df['MIT specificity']) #bigger is better
        x3 = list(pegRNA_df['CFD efficiency percentile'])

        x4 = np.array(pegRNA_df['PBS GC content']) #want this between 0.4 to 0.6

        #trnasformation that turns vals between 0.4 and 0.6 to 0
        #and penalizes significant deviations from these values
        #such that GC content = 0, or = 1 results in a value of 1 (higher = worse)
        x4_new = []
        for i in x4:
            if i<=0.6 and i>=0.4:
                x4_new.append(0)
            else:
                if i<0.4:
                    dist = 0.4-i
                    val = dist/0.4 #if GC content = 0; val = 1


                elif i>0.6:

                    dist = i-0.6
                    val = dist/0.4 #if GC content = 0; val = 1

                x4_new.append(val)


        #boolean factors
        x5 = list(pegRNA_df['contains terminator?']) #want this to be no
        x6 = list(pegRNA_df['contains RE site?']) #want this to be no
        x7 = list(pegRNA_df['last templated base C?']) #want this to be no

        x5 = np.array([int(i=='yes') for i in x5]) #turning yes/no list into 1 = yes; 0 = no
        x6 = np.array([int(i=='yes') for i in x6]) #turning yes/no list into 1 = yes; 0 = no
        x7 = np.array([int(i=='yes') for i in x7]) #turning yes/no list into 1 = yes; 0 = no

        #hard-coded weightings
        a1 = 2
        a2 = 1
        a3 = 1
        a4 = -2  #penalize deviations far from 0.4 to 0.6 GC content

        #if present, -1; otherwise, no contribution to score
        a5 = -2 #penalize doubly for terminator 
        a6 = -1
        a7 = -1



        #for distance to end of RTT from mutation, make it favored up to a distance of 15
        g = np.clip(x1/15,0,1)

        composite_score = []
        for i in range(len(g)):
            if x2[i]=='none': #exclude from scoring...could add 1 to assume about average...?
                if RE_filter==True:
                    s = a1*g[i] + a4*np.asarray(x4_new[i]) + a5*x5[i] + a6*x6[i] + a7*x7[i]
                elif RE_filter==False:
                    s = a1*g[i] + a4*np.asarray(x4_new[i]) + a5*x5[i] + a7*x7[i]



            else:
                if RE_filter==True:
                    s = a1*g[i] + a2*(x2[i]/100) + a3*(x3[i]/100) + a4*np.asarray(x4_new[i]) + a5*x5[i] + a6*x6[i] + a7*x7[i]
                elif RE_filter==False:
                    s = a1*g[i] + a2*(x2[i]/100) + a3*(x3[i]/100) + a4*np.asarray(x4_new[i]) + a5*x5[i] + a7*x7[i]

            composite_score.append(s)

        pegRNA_df['composite score']= composite_score
        
        
    
    else:


        #numerical factors
        x1 = np.array(pegRNA_df["distance mut to 5' RTT"]) #want this to be large (edit close to nick site)
        x2 = np.array(pegRNA_df['MIT specificity'])/100 #bigger is better
        x3 = np.array(pegRNA_df['CFD efficiency percentile'])/100

        x4 = np.array(pegRNA_df['PBS GC content']) #want this between 0.4 to 0.6

        #trnasformation that turns vals between 0.4 and 0.6 to 0
        #and penalizes significant deviations from these values
        #such that GC content = 0, or = 1 results in a value of 1 (higher = worse)
        x4_new = []
        for i in x4:
            if i<=0.6 and i>=0.4:
                x4_new.append(0)
            else:
                if i<0.4:
                    dist = 0.4-i
                    val = dist/0.4 #if GC content = 0; val = 1


                elif i>0.6:

                    dist = i-0.6
                    val = dist/0.4 #if GC content = 0; val = 1

                x4_new.append(val)


        #boolean factors
        x5 = list(pegRNA_df['contains terminator?']) #want this to be no
        x6 = list(pegRNA_df['contains RE site?']) #want this to be no
        x7 = list(pegRNA_df['last templated base C?']) #want this to be no

        x5 = np.array([int(i=='yes') for i in x5]) #turning yes/no list into 1 = yes; 0 = no
        x6 = np.array([int(i=='yes') for i in x6]) #turning yes/no list into 1 = yes; 0 = no
        x7 = np.array([int(i=='yes') for i in x7]) #turning yes/no list into 1 = yes; 0 = no

        #hard-coded weightings
        a1 = 2
        a2 = 1
        a3 = 1
        a4 = -2  #penalize deviations far from 0.4 to 0.6 GC content

        #if present, -1; otherwise, no contribution to score
        a5 = -2
        a6 = -1
        a7 = -1



        #for distance to end of RTT from mutation, make it favored up to a distance of 20
        g = np.clip(x1/15,0,1)
        if RE_filter==True:
            composite_score = a1*g + a2*x2 + a3*x3 + a4*np.asarray(x4_new) + a5*x5 + a6*x6 + a7*x7
        elif RE_filter==False:
            composite_score = a1*g + a2*x2 + a3*x3 + a4*np.asarray(x4_new) + a5*x5 + a7*x7

        pegRNA_df['composite score']= composite_score
    
    return pegRNA_df


def filtration(pegRNA_df, chrom_dict, RE_filter=True):
    """
    Calculates properties of pegRNAs in pegRNA_df using on_off_target(), other_filtration(), and composite_score().
    Returns pegRNA_df with these properties added.

    Parameters
    ----------
    pegRNA_df
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()
    
    chrom_dict
        *type = dict*
        
        A dictionary containing the on/off-target scores for each of the exonic regions in each chromosome.
        Keys = chromosome identifiers. See documentation for more information and access to provided reference file.
        
        If no reference file available SET THIS EQUAL TO 'none' and ranking will proceed without on/off target score.
    
        RE_filter
        *type = Boolean (True or False)

        If RE_filter = True, include Esp3I and EcoRI site inclusion in ranking system (RE sites in default adapter sequences). 
        If RE_filter= False, do not include this in the composite score...
    """
    pegRNA_df = on_off_target(pegRNA_df, chrom_dict)
    pegRNA_df = other_filtration(pegRNA_df)
    pegRNA_df = composite_score(pegRNA_df, RE_filter=RE_filter)
    
    
    
    return pegRNA_df


def ranking_system(pegRNA_df, chrom_dict, rank_by = "composite score", guides_per_mut = 5, RE_filter=True):
    
    """
    Calculates properties of pegRNAs in pegRNA_df and then ranks pegRNAs and returns the desired number of top ranked guides.
    

    Parameters
    ----------
    pegRNA_df
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()
    
    chrom_dict
        *type = dict*
        
        A dictionary containing the on/off-target scores for each of the exonic regions in each chromosome.
        Keys = chromosome identifiers. See documentation for more information and access to provided reference file.
        
        If no reference file available SET THIS EQUAL TO 'none' and ranking will proceed without on/off target score.
        
    rank_by
        *type = str*
        
        This specifies what property of the pegRNAs to use for ranking purpose (i.e. what column of pegRNA_df). This is automatically set to "composite score"
        which takes into account all of the pegRNA properties. Can alternatively be set to any of the pegRNA properties (i.e. column names) calculated
        by filtration() or any user provided pegRNA properties, as long as this is a column in pegRNA_df.
        
    guides_per_mut
        *type=int*
        
        Desired number of pegRNAs to output for each mutation. Automatically set to 5. To output all possible guides, just set
        this to a large number (e.g. 1000).

    RE_filter
        *type = Boolean (True or False)

        If RE_filter = True, include Esp3I and EcoRI site inclusion in ranking system (RE sites in default adapter sequences). 
        If RE_filter= False, do not include this in the composite score...
    
    """
    
    
    #input = pegRNA_df, ranking system, and desired # of top guides per mutation (default<=5)
    
    #how to rank the guides -- default = "composite score"
    
    #need to implement these other options:
    #other options = 'MIT specificity',  'distance mut to 5' RTT', 'CFD efficiency percentile',
       #'contains terminator?', 'contains RE site?', 'last templated base C?',
       #'PBS GC content' (all column values; may require some adjusting)
    
    
    filtered = filtration(pegRNA_df, chrom_dict, RE_filter=RE_filter)
    
    unique_idxs = np.unique(np.asarray(filtered['mutant index']))
    
    df_list = []
    
    for i in unique_idxs:
        sub_df = filtered[filtered['mutant index']==i] #consider a subset of the dataframe
        sorted_df = sub_df.sort_values(by=rank_by, ascending=False)
        top_pegRNAs = sorted_df[0:guides_per_mut]
        df_list.append(top_pegRNAs)
        
    return pd.concat(df_list)
    
    
    
#----composite run function-------

def run(mutant_input, mut_idx_list, records, index_list, minus_seqs, chrom_dict, PAM, RTT_lengths, PBS_lengths, guides_per_mut=5, RE_filter=True):
    """
    Function that (a) generates all possible pegRNAs for input mutations according to input parameters and
    (b) filters and ranks these pegRNAs according to their properties.
    
    Returns a pd.DataFrame containing the pegRNAs.
    
    Parameters
    ----------
    mutant_input
        *type = pd.DataFrame*
        
        A dataframe containing the input mutations from which selections are made to generate pegRNAs.
        See documentation for precise qualities of this dataframe
     
    mut_idx_list
        *type = list*
        
        List of indeces of mutations within mutant_input that user wants to perform function on (i.e. a subset of the mutations).        
        
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records list (reference genome). See documentation for precise requirements.
      
    minus_seqs
        *type = list*
        
        List containing reverse complement of records (reference genome). Generated by minus_seq_generator function.
    
    chrom_dict
        *type = dict*
        
        A dictionary containing the on/off-target scores for each of the exonic regions in each chromosome.
        Keys = chromosome identifiers. See documentation for more information and access to provided reference file.
        
        If no reference file available SET THIS EQUAL TO 'none' and ranking will proceed without on/off target score.
    
    PAM
        *type = str*
        
        String designating what PAM sequence to search for. Not all PAM sequences exhaustively tested.
        Normally this will be equal to 'NGG' for most existing prime editors.
        
    PBS_lengths
        *type = list*
        
        Specifies a list of the primer binding sequence (PBS) lengths of the pegRNAs being generated.
        
    RTT_lengths
        *type = list*
        
        Specifies a list of the reverse transcriptase template (RTT) lengths of the pegRNAs being generated. 
    
    guides_per_mut
        *type=int*
        
        Desired number of pegRNAs to output for each mutation. Automatically set to 5. To output all possible guides, just set
        this to a large number (e.g. 1000).

    RE_filter
        *type = Boolean (True or False)

        If RE_filter = True, include Esp3I and EcoRI site inclusion in ranking system (RE sites in default adapter sequences). 
        If RE_filter= False, do not include this in the composite score...
    """
    
    
    df_list = []

    for RTT_length in RTT_lengths:
        for PBS_length in PBS_lengths:

            marked_array_list = PAM_finder_multi_index(mutant_input, PAM, RTT_length, mut_idx_list,records, index_list)
            pegRNA_df = pegRNA_generator(mutant_input, PBS_length, RTT_length, PAM, marked_array_list, mut_idx_list,records, index_list,minus_seqs)
            df_list.append(pegRNA_df)

    concat_pegRNA_df = pd.concat(df_list)

    ranked_filtered = ranking_system(concat_pegRNA_df, chrom_dict, rank_by = "composite score", guides_per_mut = guides_per_mut, RE_filter=RE_filter)
    
    #calculating/returning consequences of edit at genomic locus...
    out = mutation_consequences(mutant_input, mut_idx_list, records, index_list)
    m = dict(zip(mut_idx_list, out))
    out_final = []
    for i, val in ranked_filtered.iterrows():
        out_final.append(m[val['mutant index']])
    
    ranked_filtered['target genomic edit'] = out_final
    
    return ranked_filtered.reset_index().drop(columns='index')


#--------viz functions------------

def align_display(pegRNA_df, records, index_list):
    
    """
    A simple visualization tool for aligning the 3' extension of the pegRNAs in pegRNA_df to the reference genome.
    Useful for spot checking pegRNA designs.

    Parameters
    ----------
    pegRNA_df
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()
        To subset (e.g. look at first 2 mutations), use a slice (e.g. pegRNA_df[0:2]).
    
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records list (reference genome). See documentation for precise requirements.
    
    """
    
    for i in range(len(pegRNA_df)):
        #-----------loading in correct mutant-----------------##
    
        #idx = pegRNA_df.iloc[[i]]['mutant index']
        mut = pegRNA_df.iloc[[i]]
        mut_type = mut['variant_type'].values[0]

                #getting gene and start/end position of mutation
        gene = mut['gene'].values[0]

        s = mut['start'].values[0]
        e = mut['end'].values[0]
        size_mut = (e-s)+1 #size of mutation; not correct for insertions though

                    ##loading in genome information
        ref_allele = mut['ref_allele'].values[0]
        mut_allele = mut['mut_allele'].values[0]

        seq_start = mut['start'].values[0]
        seq_end = mut['end'].values[0]
        chromosome = mut['chr'].values[0]
        chromosome = chromosome.split('chr')[1]

        if chromosome=='X':
            chrom = 22
        else:
            chrom = int(chromosome)-1

        seq1 = records[index_list[chrom]].seq
        
        
        
        #---------loading in information from df about pegRNA design#------------------
        protospacer = Bio.Seq.Seq(pegRNA_df.iloc[[i]]['protospacer'].values[0])
        RT = Bio.Seq.Seq(pegRNA_df.iloc[[i]]['RTT'].values[0])
        PBS = Bio.Seq.Seq(pegRNA_df.iloc[[i]]['PBS'].values[0])
        PBS_RT = Bio.Seq.Seq(pegRNA_df.iloc[[i]]['PBS_RTT_5to3'].values[0])

        ref_seq = seq1[seq_start-50:seq_start+50]

        strand = pegRNA_df.iloc[[i]]['PAM strand'].values[0]

        if strand == '+':
            alignments = pairwise2.align.localms(ref_seq.upper(), PBS_RT.reverse_complement(), 2,0, -3, -0.2)
            print(format_alignment(*alignments[0], full_sequences=True))

            #alignments = pairwise2.align.localms(ref_seq.upper(), PBS.complement(), 2,0, -3, -0.2)
            #print(format_alignment(*alignments[0], full_sequences=True))

            #alignments = pairwise2.align.localms(ref_seq.upper(), RT.complement(), 2,0, -3, -0.2)
            #print(format_alignment(*alignments[0], full_sequences=True))

        else:

            alignments = pairwise2.align.localms(ref_seq.upper(), PBS_RT, 2,0, -3, -0.2)
            print(format_alignment(*alignments[0], full_sequences=True))



            #alignments = pairwise2.align.localms(ref_seq.upper(), PBS, 2,0, -3, -0.2)
            #print(format_alignment(*alignments[0], full_sequences=True))

            #alignments = pairwise2.align.localms(ref_seq.upper(), RT, 2,0, -3, -0.2)
            #print(format_alignment(*alignments[0], full_sequences=True))

            #alignments = pairwise2.align.localms(ref_seq.upper().reverse_complement(), protospacer, 2,0, -3, -0.2)
            #print(format_alignment(*alignments[0], full_sequences=True))
            
            
def split(word):
    """Simple function for splitting string into component characters"""
    return [char for char in word]


def pegrna_display(pegRNA_df, pegRNA_df_loc, records, index_list):
    """
    A visualization tool for looking at pegRNAs.
    Returns a figure that can be saved or shown inline.

    Parameters
    ----------
    pegRNA_df
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()
    
    pegRNA_df_loc
        *type = int*
        
        Which iloc (row index) in pegRNA_df_loc to generate the visualizaiton for. 
    
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records list (reference genome). See documentation for precise requirements.
    
    """
    
    i = pegRNA_df_loc
    
    #idx = pegRNA_df.iloc[[i]]['mutant index']
    mut = pegRNA_df.iloc[[i]]
    mut_type = mut['variant_type'].values[0]

            #getting gene and start/end position of mutation
    gene = mut['gene'].values[0]

    s = mut['start'].values[0]
    e = mut['end'].values[0]
    size_mut = (e-s)+1 #size of mutation; not correct for insertions though

                ##loading in genome information
    ref_allele = mut['ref_allele'].values[0]
    mut_allele = mut['mut_allele'].values[0]

    seq_start = mut['start'].values[0]
    seq_end = mut['end'].values[0]
    chromosome = mut['chr'].values[0]
    chromosome = chromosome.split('chr')[1]
    
    if chromosome=='X':
        chrom = 22
    else:
        chrom = int(chromosome)-1
    
    seq1 = records[index_list[chrom]].seq

    #---------loading in information from df about pegRNA design#------------------
    protospacer = Bio.Seq.Seq(pegRNA_df.iloc[[i]]['protospacer'].values[0])
    RT = Bio.Seq.Seq(pegRNA_df.iloc[[i]]['RTT'].values[0])
    PBS = Bio.Seq.Seq(pegRNA_df.iloc[[i]]['PBS'].values[0])
    var_type = pegRNA_df.iloc[[i]]['variant_type'].values[0]
    nick_dist = pegRNA_df.iloc[[i]]['distance to nick'].values[0]
    
    PBS_length = len(PBS)
    RTT_length = len(RT)
    
    PBS_RT = Bio.Seq.Seq(pegRNA_df.iloc[[i]]['PBS_RTT_5to3'].values[0])

    ref_seq = seq1[seq_start-50:seq_start+50].upper()

    strand = pegRNA_df.iloc[[i]]['PAM strand'].values[0]
    
    #-------display-----------------#
    dict_bases = {'T':0, 'A':1, 'C':2, 'G':3, '-':4, 'X':5}
    
    target_seq = pegRNA_df.iloc[[i]]['target sequence'].values[0]
    
    
    #for plus strand
    if strand=='+':
    
        #locating pam
        PAM_n = pegRNA_df.iloc[[i]]['PAM location'].values[0] - pegRNA_df.iloc[[i]]['start'].values[0] + 50

        
        #protospacer processing
        proto = pegRNA_df.iloc[[i]]['protospacer'].values[0]
        protospacer_l = ['-']*(PAM_n - 20)
        protospacer_r = ['-']* (100-PAM_n)
        proto_middle = split(proto)
        #print(len(protospacer_l) + len(protospacer_r) + len(proto_middle ))
        protospacer = protospacer_l+proto_middle+protospacer_r

        #target_sequence_processing
        target_seq_l = ['-']*(PAM_n - 20-5)
        target_seq_r = ['-']* (100-(PAM_n + (len(target_seq)-20))+5)
        target_seq_mid = split(target_seq)
        target = target_seq_l + target_seq_mid + target_seq_r
        
        ##and the reference sequence
        test = str(ref_seq)
        test_comp = Bio.Seq.Seq(test).complement()

        split_test = split(test)
        split_test_comp = split(str(test_comp))

        num_translation = [dict_bases[i] for i in split_test]
        num_translation_c = [dict_bases[i] for i in split_test_comp]
        
        
        #pbs rtt processing; display varies according to type...
        misc = ['SNP','ONP','DNP']
        if var_type in misc:
        
            pbs_rtt_5to3 = pegRNA_df.iloc[[i]]['PBS_RTT_5to3'].values[0]
            pbs_rtt_3to5 = str(Bio.Seq.Seq(pbs_rtt_5to3).reverse_complement().complement())

            RTT_length = len(pbs_rtt_3to5)-PBS_length
            pbs_rtt_start = PAM_n-3-PBS_length

            PBS_RTT = split(pbs_rtt_3to5)
            PBS_RTT_l = ['-']*(pbs_rtt_start)
            PBS_RTT_r = ['-']*(100-(RTT_length+PBS_length+pbs_rtt_start))
            PBS_RTT_true = PBS_RTT_l + PBS_RTT + PBS_RTT_r
            
        elif var_type=='DEL':
            
            if mut_allele == '-':
                diff = len(ref_allele)-len(mut_allele)+1
                
                pbs_rtt_5to3 = pegRNA_df.iloc[[i]]['PBS_RTT_5to3'].values[0]
                pbs_rtt_3to5 = str(Bio.Seq.Seq(pbs_rtt_5to3).reverse_complement().complement())

                RTT_length = len(pbs_rtt_3to5)-PBS_length
                pbs_rtt_start = PAM_n-3-PBS_length

                PBS_RTT = split(pbs_rtt_3to5)
                PBS_RTT_l = ['-']*(pbs_rtt_start)
                PBS_RTT_r = ['-']*(100-(RTT_length+PBS_length+pbs_rtt_start+diff))

                PBS_RTT_true = PBS_RTT_l + PBS_RTT[0:PBS_length+nick_dist-1] + ['X']*diff + PBS_RTT[PBS_length+nick_dist-1:]+ PBS_RTT_r

                
            else:
                diff = len(ref_allele)-len(mut_allele)
            
                pbs_rtt_5to3 = pegRNA_df.iloc[[i]]['PBS_RTT_5to3'].values[0]
                pbs_rtt_3to5 = str(Bio.Seq.Seq(pbs_rtt_5to3).reverse_complement().complement())

                RTT_length = len(pbs_rtt_3to5)-PBS_length
                pbs_rtt_start = PAM_n-3-PBS_length

                PBS_RTT = split(pbs_rtt_3to5)
                PBS_RTT_l = ['-']*(pbs_rtt_start)
                PBS_RTT_r = ['-']*(100-(RTT_length+PBS_length+pbs_rtt_start+diff))

                PBS_RTT_true = PBS_RTT_l + PBS_RTT[0:PBS_length+nick_dist] + ['X']*diff + PBS_RTT[PBS_length+nick_dist:]+ PBS_RTT_r
            
        elif var_type=='INS':
            
            pbs_rtt_5to3 = pegRNA_df.iloc[[i]]['PBS_RTT_5to3'].values[0]
            pbs_rtt_3to5 = str(Bio.Seq.Seq(pbs_rtt_5to3).reverse_complement().complement())

            RTT_length = len(pbs_rtt_3to5)-PBS_length
            pbs_rtt_start = PAM_n-3-PBS_length

            PBS_RTT = split(pbs_rtt_3to5)
            PBS_RTT_l = ['-']*(pbs_rtt_start)
            PBS_RTT_r = ['-']*(100-(RTT_length+PBS_length+pbs_rtt_start))
            PBS_RTT_true = PBS_RTT_l + PBS_RTT + PBS_RTT_r


        
        blank1 = ['-']*100
        blank1_trans = [dict_bases[i] for i in blank1]
        
        text_df = [PBS_RTT_true, split_test, split_test_comp, protospacer, blank1, target]

        proto_translation = [dict_bases[i] for i in protospacer]
        PBS_rtt_trans = [dict_bases[i] for i in PBS_RTT_true]
        target_translation = [dict_bases[i] for i in target]

        dataFrame = [PBS_rtt_trans,num_translation,num_translation_c, proto_translation, blank1_trans, target_translation]

        #-------plotting-------------------#
        # For only three colors, it's easier to choose them yourself.
        # If you still really want to generate a colormap with cubehelix_palette instead,
        # add a cbar_kws={"boundaries": linspace(-1, 1, 4)} to the heatmap invocation
        # to have it generate a discrete colorbar instead of a continous one.
        fig_height = 4.5
        fig = plt.figure(figsize=(20,fig_height))
        myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white')
        if var_type=='DEL':
            myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white', 'black')

        cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

        ax = sns.heatmap(dataFrame, annot=text_df, fmt="", linewidth=0, cmap=cmap, cbar=False,linewidths=.5, linecolor='lightgray')

        # Manually specify colorbar labelling after it's been generated
        #colorbar = ax.collections[0].colorbar
        #height = len(myColors)
        #colorbar.set_ticks([height/5 - height/8 , 2*height/5 - height/8, 
        #                    3*height/5 - height/8, 4*height/5 - height/8,5*height/5 - height/8,  ])

        #colorbar.set_ticks([0,1,2,3,4])
        #colorbar.set_ticklabels(['T', 'A','C', 'G', '-'])

        # X - Y axis labels
        #ax.set_ylabel('FROM')
        #ax.set_xlabel('TO')

        #add patches
        #PAM patch
        ax.add_patch(patches.Rectangle((PAM_n, 1), 3, 1, fill=False, edgecolor='tab:green', lw=3, label='PAM')) #protospacer loca)) #protospacer location

        #target edit patch
        s = pegRNA_df.iloc[[i]]['start'].values[0]
        e = pegRNA_df.iloc[[i]]['end'].values[0]
        ax.add_patch(patches.Rectangle((49, 1), e-s+1, 1, fill=False, edgecolor='tab:red', lw=3, label='target edit site')) #protospacer location

        #protospacer
        ax.add_patch(patches.Rectangle((PAM_n-20, 3), 20, 1, fill=False, edgecolor='tab:blue', lw=3, label='protospacer')) #protospacer location

        #protospacer
        ax.add_patch(patches.Rectangle((pbs_rtt_start, 0), PBS_length, 1, fill=False, edgecolor='tab:purple', lw=3, label='PBS')) #protospacer location

        #RTT
        
        
        if var_type=='DEL':
            ax.add_patch(patches.Rectangle((pbs_rtt_start+PBS_length, 0), RTT_length+diff, 1, fill=False, edgecolor='yellow', lw=3, label='RTT')) #protospacer location
            ax.add_patch(Polygon([(pbs_rtt_start+PBS_length+RTT_length+diff, 0), (pbs_rtt_start+PBS_length+RTT_length+diff+2, 0.5), (pbs_rtt_start+PBS_length+RTT_length+diff, 1)],facecolor='black',edgecolor='yellow', lw=3))
        #ax.add_patch(Polygon([(70, 0), (72, 0.5), (70, 1)],edgecolor='yellow',fill='yellow', lw=3))


        
        else:
            ax.add_patch(patches.Rectangle((pbs_rtt_start+PBS_length, 0), RTT_length, 1, fill=False, edgecolor='yellow', lw=3, label='RTT')) #protospacer location
            ax.add_patch(Polygon([(pbs_rtt_start+PBS_length+RTT_length, 0), (pbs_rtt_start+PBS_length+RTT_length+2, 0.5), (pbs_rtt_start+PBS_length+RTT_length, 1)],facecolor='black',edgecolor='yellow', lw=3))



        #target sequence
        ax.add_patch(patches.Rectangle((PAM_n - 25, 5), len(target_seq), 1, fill=False, edgecolor='black', lw=3, label='Target seq')) #protospacer location

        
        ax.set_title('RTT length: '+str(RTT_length)+ ' nt, PBS length: ' + str(PBS_length)+' nt | ' + mut_type+ ': '+ref_allele + '>' + mut_allele, fontsize=16)


        ax.set_yticklabels(['3', '5', '3', '5','', '5'])

        # Only y-axis labels need their rotation set, x-axis labels already have a rotation of 0
        _, labels = plt.yticks()
        plt.setp(labels, rotation=0)

        plt.legend(bbox_to_anchor=(1., 1.), loc='upper left', fontsize=20)
        plt.tight_layout()
        #plt.show()

    if strand=='-':
    
        #locating pam
        PAM_n = pegRNA_df.iloc[[i]]['PAM location'].values[0] - pegRNA_df.iloc[[i]]['start'].values[0]+50-1

        #protospacer processing
        proto = pegRNA_df.iloc[[i]]['protospacer'].values[0]
        proto = str(Bio.Seq.Seq(proto).reverse_complement().complement())
        protospacer_l = ['-']*(PAM_n+1)
        protospacer_r = ['-']* (100-PAM_n-21)
        proto_middle = split(proto)
        #print(len(protospacer_l) + len(protospacer_r) + len(proto_middle ))
        protospacer = protospacer_l+proto_middle+protospacer_r

        #target_sequence_processing
        target_seq_l = ['-']*(PAM_n+1-(len(target_seq)-20)+5)
        target_seq_r = ['-']* (100-21-(PAM_n)-5)
        target_seq_mid = split(target_seq)
        target = target_seq_l + target_seq_mid + target_seq_r
        
        
        #pbs rtt processing
                #pbs rtt processing; display varies according to type...
        misc = ['SNP','ONP','DNP']
        if var_type in misc:
        
            pbs_rtt_5to3 = pegRNA_df.iloc[[i]]['PBS_RTT_5to3'].values[0]
            #pbs_rtt_3to5 = str(Bio.Seq.Seq(pbs_rtt_5to3).reverse_complement().complement())

            RTT_length = len(pbs_rtt_5to3)-PBS_length
            pbs_rtt_start = PAM_n+4+PBS_length

            PBS_RTT = split(pbs_rtt_5to3)
            PBS_RTT_r = ['-']*(100-pbs_rtt_start)
            PBS_RTT_l = ['-']*(pbs_rtt_start- (RTT_length+PBS_length))
            PBS_RTT_true = PBS_RTT_l + PBS_RTT + PBS_RTT_r
            
        elif var_type=='DEL':
            
            if mut_allele == '-':
                diff = len(ref_allele)-len(mut_allele)+1
                
            else:
                diff = len(ref_allele)-len(mut_allele)
            
            
            pbs_rtt_5to3 = pegRNA_df.iloc[[i]]['PBS_RTT_5to3'].values[0]
            #pbs_rtt_3to5 = str(Bio.Seq.Seq(pbs_rtt_5to3).reverse_complement().complement())

            RTT_length = len(pbs_rtt_5to3)-PBS_length
            pbs_rtt_start = PAM_n+4+PBS_length

            PBS_RTT = split(pbs_rtt_5to3)
            PBS_RTT_r = ['-']*(100-pbs_rtt_start)
            PBS_RTT_l = ['-']*(pbs_rtt_start- (RTT_length+PBS_length+diff))
            PBS_RTT_true = PBS_RTT_l + PBS_RTT + PBS_RTT_r
            
            PBS_RTT_true = PBS_RTT_l + PBS_RTT[0:(RTT_length-nick_dist)] + ['X']*diff + PBS_RTT[(RTT_length-nick_dist):]+ PBS_RTT_r
            
        elif var_type=='INS':
            
            pbs_rtt_5to3 = pegRNA_df.iloc[[i]]['PBS_RTT_5to3'].values[0]
            pbs_rtt_3to5 = str(Bio.Seq.Seq(pbs_rtt_5to3).reverse_complement().complement())

            RTT_length = len(pbs_rtt_3to5)-PBS_length
            pbs_rtt_start = PAM_n-3-PBS_length

            PBS_RTT = split(pbs_rtt_3to5)
            PBS_RTT_l = ['-']*(pbs_rtt_start)
            PBS_RTT_r = ['-']*(100-(RTT_length+PBS_length+pbs_rtt_start))
            PBS_RTT_true = PBS_RTT_l + PBS_RTT + PBS_RTT_r
        
        
        ##and the reference sequence
        test = str(ref_seq)
        test_comp = Bio.Seq.Seq(test).complement()

        split_test = split(test)
        split_test_comp = split(str(test_comp))

        num_translation = [dict_bases[i] for i in split_test]
        num_translation_c = [dict_bases[i] for i in split_test_comp]

        blank1 = ['-']*100
        blank1_trans = [dict_bases[i] for i in blank1]
        
        text_df = [ protospacer, split_test, split_test_comp,PBS_RTT_true,blank1, target]

        target_translation = [dict_bases[i] for i in target]
        proto_translation = [dict_bases[i] for i in protospacer]
        PBS_rtt_trans = [dict_bases[i] for i in PBS_RTT_true]

        dataFrame = [proto_translation, num_translation,num_translation_c,  PBS_rtt_trans,blank1_trans,target_translation]

        #-------plotting-------------------#
        # For only three colors, it's easier to choose them yourself.
        # If you still really want to generate a colormap with cubehelix_palette instead,
        # add a cbar_kws={"boundaries": linspace(-1, 1, 4)} to the heatmap invocation
        # to have it generate a discrete colorbar instead of a continous one.
        fig_height = 4.5
        fig = plt.figure(figsize=(20,fig_height))
        myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white')
        if var_type=='DEL':
            myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white', 'black')

        cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

        ax = sns.heatmap(dataFrame, annot=text_df, linewidth=0, fmt="", cmap=cmap, cbar=False, linewidths=.5, linecolor='lightgray')

        # Manually specify colorbar labelling after it's been generated

        # X - Y axis labels
        #ax.set_ylabel('FROM')
        ax.set_yticklabels(['3', '5', '3', '5','', '5'])

        #add patches
        #PAM patch
        ax.add_patch(patches.Rectangle((PAM_n-2, 2), 3, 1, fill=False, edgecolor='tab:green', lw=3, label='PAM')) #protospacer location

        #target edit patch
        s = pegRNA_df.iloc[[i]]['start'].values[0]
        e = pegRNA_df.iloc[[i]]['end'].values[0]
        ax.add_patch(patches.Rectangle((49, 1), e-s+1, 1, fill=False, edgecolor='tab:red', lw=3, label='target edit site')) #protospacer location

        #protospacer
        ax.add_patch(patches.Rectangle((PAM_n+1, 0), 20, 1, fill=False, edgecolor='tab:blue', lw=3, label='protospacer')) #protospacer location

        #PBS
        ax.add_patch(patches.Rectangle((pbs_rtt_start-PBS_length, 3), PBS_length, 1, fill=False, edgecolor='tab:purple', lw=3, label='PBS')) #protospacer location

        #RTT
        if var_type=='DEL':
            #ax.add_patch(patches.Rectangle((pbs_rtt_start+PBS_length, 0), RTT_length+diff, 1, fill=False, edgecolor='yellow', lw=3, label='RTT')) #protospacer location
            
            ax.add_patch(patches.Rectangle((pbs_rtt_start-PBS_length-RTT_length-diff, 3), RTT_length+diff, 1, fill=False, edgecolor='yellow', lw=3, label='RTT')) #protospacer location
            
            ax.add_patch(Polygon([(pbs_rtt_start-PBS_length-RTT_length-diff, 3), (pbs_rtt_start-PBS_length-RTT_length-diff-2, 3.5), (pbs_rtt_start-PBS_length-RTT_length-diff, 4)],fill='yellow',facecolor='black',edgecolor='yellow', lw=3))



        else:
            ax.add_patch(patches.Rectangle((pbs_rtt_start-PBS_length-RTT_length, 3), RTT_length, 1, fill=False, edgecolor='yellow', lw=3, label='RTT')) #protospacer location
            ax.add_patch(Polygon([(pbs_rtt_start-PBS_length-RTT_length, 3), (pbs_rtt_start-PBS_length-RTT_length-2, 3.5), (pbs_rtt_start-PBS_length-RTT_length, 4)],facecolor='black',edgecolor='yellow', lw=3))




        
        
        #target sequence
        ax.add_patch(patches.Rectangle(((PAM_n+1-(len(target_seq)-20))+5, 5), len(target_seq), 1, fill=False, edgecolor='black', lw=3, label='Target seq')) #protospacer location

        
        ax.set_title('RTT length: '+str(RTT_length)+ ' nt, PBS length: ' + str(PBS_length)+' nt | ' + mut_type+ ': '+ref_allele + '>' + mut_allele, fontsize=16)

        plt.legend(bbox_to_anchor=(1., 1.), loc='upper left', fontsize=20)
        # Only y-axis labels need their rotation set, x-axis labels already have a rotation of 0
        _, labels = plt.yticks()
        plt.setp(labels, rotation=0)
        plt.tight_layout()
        #plt.show()
        
    return fig

#----------automated oligo generation function---------

def oligo_generator(ranked_filtered,five_prime_adapter = 'AGCGTACACGTCTCACACC',three_prime_adapter = 'GAATTCTAGATCCGGTCGTCAAC',
                    gRNA_scaff = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC',
                    sensor=True, append_proto_G = True):
    #generates oligo sequences based on generate pegRNA sequences
    #if sensor=True, append synthetic version of target sequence within oligo
    #if sensor=False, only generate pegRNA (no synthetic target sequence)
    """
    A tool for automatically generating oligos from the output of run().
    Returns input dataframe with new columns containing pegRNA version and epegRNA version of oligo.
    (epegRNA version just contains additional 3' structural motif, tevopreQ1, which helps prevent degradation
    of the pegRNA)

    Parameters
    ----------
    ranked_filtered
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()
    
    five_prime_adapter
        *type = str*
        
        5' Prime Adapter. The automatically provided 5' adapter contains an Esp3I (BsmBI) site. Can be swapped with 
        whatever input string user wants.
    
    three_prime_adapter
        *type = str*
        
        5' Prime Adapter. The automatically provided 5' adapter contains an Esp3I (BsmBI) site. Can be swapped with 
        whatever input string user wants.
        
    gRNA_scaff
        *type = str*
        
        gRNA scaffold region. Automatically set to a functional gRNA scaffold. Can be swapped with 
        whatever input string user wants.
        
    sensor
        *type = boolean (True or False)*
        
        If true, include synthetic version of endogenous target site ("sensor" region) in oligo.
        See documentation for more details.
        
    append_proto_G
        *type = boolean (True or False)*
        
        If true, append a 'G' to the beginning of the protospacer. Reccomended for efficient transcription of pegRNA.
        
    
    """
        
        
    appended_proto_G = 'G'
    u6_term = 'TTTTTTT'
    tevopreQ1 = 'CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA'
    
    protos = list(ranked_filtered['protospacer'])
    pbs_rtt = list(ranked_filtered['PBS_RTT_5to3'])
    target_seq = list(ranked_filtered['target sequence'])
    
    peg = []
    epeg = []
        
    if append_proto_G == True:
        if sensor==True:
            for i in range(len(ranked_filtered)):

                oli = five_prime_adapter+appended_proto_G+protos[i]+gRNA_scaff+pbs_rtt[i]+u6_term+target_seq[i]+three_prime_adapter

                peg.append(oli)

                oli_e = five_prime_adapter+appended_proto_G+protos[i]+gRNA_scaff+pbs_rtt[i]+tevopreQ1+u6_term+target_seq[i]+three_prime_adapter

                epeg.append(oli_e)
                
        elif sensor==False:
            for i in range(len(ranked_filtered)):

                oli = five_prime_adapter+appended_proto_G+protos[i]+gRNA_scaff+pbs_rtt[i]+u6_term+three_prime_adapter

                peg.append(oli)

                oli_e = five_prime_adapter+appended_proto_G+protos[i]+gRNA_scaff+pbs_rtt[i]+tevopreQ1+u6_term+three_prime_adapter

                epeg.append(oli_e)
            
    
    elif append_proto_G==False:
        if sensor==True:
            for i in range(len(ranked_filtered)):

                oli = five_prime_adapter+protos[i]+gRNA_scaff+pbs_rtt[i]+u6_term+target_seq[i]+three_prime_adapter

                peg.append(oli)

                oli_e = five_prime_adapter+protos[i]+gRNA_scaff+pbs_rtt[i]+tevopreQ1+u6_term+target_seq[i]+three_prime_adapter

                epeg.append(oli_e)
        elif sensor==False:
            for i in range(len(ranked_filtered)):

                oli = five_prime_adapter+protos[i]+gRNA_scaff+pbs_rtt[i]+u6_term+three_prime_adapter

                peg.append(oli)

                oli_e = five_prime_adapter+protos[i]+gRNA_scaff+pbs_rtt[i]+tevopreQ1+u6_term+three_prime_adapter

                epeg.append(oli_e)
    
    ranked_filtered['pegRNA_oligo']=peg
    ranked_filtered['epegRNA_tevopreQ1_oligo']=epeg
        
    
    return ranked_filtered


#----------automated library generation functions----------------
def neutral_substitutions(gene_name, chrom, strand, start_end_cds, records, index_list):
    """
    A function for generating all possible synonymous codon substitutions (i.e. silent mutations) of a given gene.
    See documentation for more information about format of start_end_cds & example usage.

    Parameters
    ----------
    gene_name
        *type = str*
        
        Gene's Hugo_Symbol (i.e. name).
        
    chrom
        *type = str*
        
        Chromosome that gene occurs on. Format = e.g. 'chr17', 'chrX', etc.
        
    strand
        *type = str*
        
        Strand that transcript is on. Options are '+' or '-'.
        
    start_end_cds
        *type = list*
        
        A 2-d list containing the start/end locations of each region of the coding sequence (CDS) for the gene's selected transcript.
        See documentation for example and precise specifications of format. 
        
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records list (reference genome). See documentation for precise requirements.
    
    
    """
    
    
    codons = []
    for i in start_end_cds:
        for k in range(i[0], i[1]+1):
            codons.append(k)

    gene_codons = []
    if strand=='+':
        for i in range(0,len(codons),3):
            gene_codons.append([codons[i], codons[i+1], codons[i+2]])

    elif strand=='-':
        codons = codons[::-1]
        for i in range(0,len(codons),3):
            gene_codons.append([codons[i], codons[i+1], codons[i+2]])
    
    
    #-------------------creating df of codons  
    if chrom=='chrX':
        chrom=23
    else:
        chrom = int(re.findall(r'\d+', chrom)[0])
    
    #loading in reference genome from peg engine module
    seq1 = records[index_list[chrom-1]].seq

    codons = []
    for i in gene_codons:
        s1 = seq1[i[0]-1]+seq1[i[1]-1]+seq1[i[2]-1]
        se = Bio.Seq.Seq(s1)
        codons.append(str(se))



    codon_start = []
    codon_end = []
    ref_codon = []
    hugo_symbol = [gene_name]*len(gene_codons)
    chrom1 = [chrom]*len(gene_codons)
    variant_type = ['ONP']*len(gene_codons)
    strand1 = [strand]*len(codons)

    for idx, val in enumerate(gene_codons):

        if strand=='-':
            ref_codon.append(codons[idx][::-1]) #need it in + strand orientation

        else:
            ref_codon.append(codons[idx]) #need it in + strand orientation


        s = val[2] #start = end
        e = val[0]
        codon_start.append(s)
        codon_end.append(e)


    g_cod = pd.DataFrame(data = hugo_symbol, columns = ['Hugo_Symbol'])
    g_cod['Chromosome']=chrom1
    g_cod['Start_Position'] = codon_start
    g_cod['End_Position']=codon_end
    g_cod['Variant_Type'] = variant_type
    g_cod['Reference_Allele'] = ref_codon
    g_cod['Tumor_Seq_Allele2']=ref_codon #change to mutations later
    g_cod['codon']=range(1,len(gene_codons)+1)


    #it's a minus strand protein
    if strand=='-':

        codon_list = list(g_cod['Reference_Allele'])
        codon_list_rc = [str(Bio.Seq.translate(Bio.Seq.transcribe(Bio.Seq.Seq(i).reverse_complement()))) for i in codon_list]

    else:

        codon_list = list(g_cod['Reference_Allele'])
        codon_list_rc = [str(Bio.Seq.translate(Bio.Seq.transcribe(Bio.Seq.Seq(i)))) for i in codon_list]


    g_cod['ref_aa']=codon_list_rc
    
    
    #------------------enumerate all synonymous codon substitutions at each position
    all_codons = []
    bases = ['A','T','C','G']
    for i in bases:
        for k in bases:
            for j in bases:
                new = i+k+j
                all_codons.append(new)


    #translations for all codons
    aas=[str(Bio.Seq.Seq(i).transcribe().translate()) for i in all_codons]
    unique_aas = list(np.unique(aas))

    #creating list corresponding to codons for each amino acid
    matched_codons = [[] for x in range(len(unique_aas))]
    index_list = []
    for idx, val in enumerate(all_codons):
        aa = str(Bio.Seq.Seq(val).transcribe().translate()) #determine amino acid from codon
        index=unique_aas.index(aa) #find index of aa in unique_aa list
        matched_codons[index].append(val)



    #------------------now creating one synonymous subtitution for each aa (iff possible)
    tumor_allele = []
    mut_aa = []
    for idx, val in enumerate(list(g_cod['ref_aa'])):

        ref_cod = g_cod.iloc[[idx]]['Reference_Allele'].values[0]


        #finding corresponding codons
        cod_idx = unique_aas.index(val)
        possible_codons = matched_codons[cod_idx]

        if strand=='-':

            possible_codons = [str(Bio.Seq.Seq(i).reverse_complement()) for i in possible_codons]

            if len(possible_codons)==1:
                tumor_allele.append(possible_codons[0])
                mut_aa.append(str(Bio.Seq.Seq(possible_codons[0]).reverse_complement().transcribe().translate()))

            else:
                if possible_codons[0]==ref_cod:
                    tumor_allele.append(possible_codons[1])
                    mut_aa.append(str(Bio.Seq.Seq(possible_codons[1]).reverse_complement().transcribe().translate()))



                else:
                    tumor_allele.append(possible_codons[0])
                    mut_aa.append(str(Bio.Seq.Seq(possible_codons[0]).reverse_complement().transcribe().translate()))


        else:

            if len(possible_codons)==1:
                tumor_allele.append(possible_codons[0])
                mut_aa.append(str(Bio.Seq.Seq(possible_codons[0]).transcribe().translate()))


            else:

                if possible_codons[0]==ref_cod:
                    tumor_allele.append(possible_codons[1])
                    mut_aa.append(str(Bio.Seq.Seq(possible_codons[1]).transcribe().translate()))


                else:
                    tumor_allele.append(possible_codons[0])
                    mut_aa.append(str(Bio.Seq.Seq(possible_codons[0]).transcribe().translate()))





    g_cod['Tumor_Seq_Allele2']=tumor_allele
    g_cod['mut_aa']=mut_aa

    #------------------classifying and filtering

    #now classifying mutations according to original and mutant
    classifier = []

    for i, val in g_cod.iterrows():

        o = val['ref_aa']
        m = val['mut_aa']
        o_dna = val['Reference_Allele']
        m_dna = val['Tumor_Seq_Allele2']

        if o_dna == m_dna:
            classifier.append('unchanged')

        else:
            if o==m:
                classifier.append('neutral')
            elif o != m:
                if m=='*':
                    classifier.append('stop')
                else:
                    classifier.append('missense')

    g_cod['classification']=classifier
    return g_cod[g_cod['classification']=='neutral'].reset_index().drop(columns='index')


def mutation_aggregator(mutant_input, gene_name):    
    """
    Selects the mutations that correspond to the desired gene name.
    Removes duplicates (checking at DNA level; not amino acid level).
    Returns dataframe containing these aggregated mutants occuring in gene_name.
    
    Parameters
    ----------
    mutant_input
        *type = pd.DataFrame*
        
        A dataframe containing the input mutations from which selections are made to generate pegRNAs.
        See documentation for precise qualities of this dataframe
        
    gene_name
        *type = str*
        
        Gene's Hugo_Symbol (i.e. name).
        
    """
    
    gene_mutants = mutant_input[mutant_input['Hugo_Symbol']==gene_name]
    mutant_input_sparse = gene_mutants[["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2"]].drop_duplicates()
    mutant_idxs = list(mutant_input_sparse.index)
    
    return mutant_input.iloc[mutant_idxs]

def library_input_generator(mutant_input, gene_name, chrom, strand, start_end_cds, records, index_list, control_fraction):
  
    """
    A function for generating a library of input mutations corresponding to all of the mutations occuring in the 
    desired gene to profile. In other words, this a tool for generating a library of input mutations for partial or
    full saturation mutagenesis of a target gene. Includes silent substitution control pegRNAs, to act as internal controls.

    Parameters
    ----------
    mutant_input
        *type = pd.DataFrame*
        
        A dataframe containing the input mutations from which selections are made to generate pegRNAs.
        See documentation for precise qualities of this dataframe
        
    gene_name
        *type = str*
        
        Gene's Hugo_Symbol (i.e. name).
        
    chrom
        *type = str*
        
        Chromosome that gene occurs on. Format = e.g. 'chr17', 'chrX', etc.
        
    strand
        *type = str*
        
        Strand that transcript is on. Options are '+' or '-'.
        
    start_end_cds
        *type = list*
        
        A 2-d list containing the start/end locations of each region of the coding sequence (CDS) for the gene's selected transcript.
        See documentation for example and precise specifications of format. 
        
    records
        *type = list*
        
        List containing reference genome. See documentation for precise requirements.
        
    index_list
        *type = list*
        
        Index for records list (reference genome). See documentation for precise requirements.
    
    control_fraction
        *type = int*
        
        Desired fraction of library that the user wants to be silent substitutions/control guides. Can range from 0 to 1.
        Reccomended ~.05 (5%).

    """
    
        
    df = mutation_aggregator(mutant_input, gene_name)
    
    #neutral_muts= neutral_substitutions(gene_name, codon_locs)
    neutral_muts = neutral_substitutions(gene_name, chrom, strand, start_end_cds, records, index_list)
    
    cod = list(neutral_muts['codon'])
    
    num_desired1 = int(len(df)*control_fraction/(1-control_fraction)) #getting number of muts for desired frac
    num_desired = min(num_desired1, len(cod))
    
    selected_cods = []
    for i in range(0,len(cod),len(cod)//num_desired):
        #iterate through the 
        selected_cods.append(cod[i])
    
    
    #subset based on selected codons
    neutral_selected = neutral_muts[neutral_muts['codon'].isin(selected_cods)]


    
    compiled_df = pd.concat((neutral_selected,df))
    
    return compiled_df.reset_index().drop(columns=['index'])


#---------library viz functions------------#
def lollipop_library(rf, gene_name, start_end_cds, strand, plot=True):
    """
    Visualization function for looking at the number of pegRNAs generated in each codon of target gene.
    Stratifies mutations by type and returns a 2-d array containing # of pegRNAs designed in each codon
    for each variant type (SNP, INS, DEL, ONP/neutral).
    
    Parameters
    ----------
    rf
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs generated by run().
        
    gene_name
        *type = str*
        
        Gene's Hugo_Symbol (i.e. name).
        
    start_end_cds
        *type = list*
        
        A 2-d list containing the start/end locations of each region of the coding sequence (CDS) for the gene's selected transcript.
        See documentation for example and precise specifications of format. 
        
    strand
        *type = str*
        
        Strand that transcript is on. Options are '+' or '-'.
        
    plot
        *type = boolean ('True' or 'False')*
        
        If true, show stacked bar plot. If false, simply return the zero_array.
        

    """
    df = rf
    
    codons = []
    for i in start_end_cds:
        for k in range(i[0], i[1]+1):
            codons.append(k)

    gene_codons = []
    if strand=='+':
        for i in range(0,len(codons),3):
            gene_codons.append([codons[i], codons[i+1], codons[i+2]])

    elif strand=='-':
        codons = codons[::-1]
        for i in range(0,len(codons),3):
            gene_codons.append([codons[i], codons[i+1], codons[i+2]])
    
    gene_codons = np.asarray(gene_codons)
    
    len_gene = len(gene_codons)
    
    zero_array = np.zeros((4,len_gene))
    
    dict_1 = {"SNP":0, "INS":1, "DEL":2,  "ONP":3}
    
    
    #iterate through mutations and determine where it occurs in gene
    
    for i, val in df.iterrows():
        
        #HGVSp = val['HGVSp']
        #cod = val['codon']
        #ref = val['ref_aa']
        #mut=val['mut_aa']
        #num_occur = val['num_occurences']
        i,j = np.where(gene_codons == val['start'])
        if len(i)>0:
            cod = i[0]
        
            var = val['variant_type']
            idx = dict_1[var]
            zero_array[idx][int(cod)]+= 1
            
        else:
            continue
       
    if plot==True:
        #plotting
        plt.figure(figsize=(12,6))
        plt.bar(range(1, len_gene+1),zero_array[0], color='tab:blue',alpha=0.7, label='SNP')
        plt.bar(range(1, len_gene+1),zero_array[1], bottom=zero_array[0], color='tab:green',alpha=0.7, label='INS')
        plt.bar(range(1, len_gene+1),zero_array[2], bottom=zero_array[0]+zero_array[1], color='tab:purple', alpha=0.7,label='DEL')
        plt.bar(range(1, len_gene+1),zero_array[3],bottom=zero_array[0]+zero_array[1]+zero_array[2], color='grey',alpha=0.7, label='Neutral')
        plt.xlim(1,len_gene+1)
        plt.legend(fontsize=15)
        plt.ylabel('Number of pegRNAs', fontsize=15)
        plt.xlabel('Codon', fontsize=15)

    
    return zero_array

def matrix_rep_library(rf,gene_name,start_end_cds, strand, plot=True):
    #requires adding HGVSp information to SNPs in dataset...
    """
    Visualization function for looking at the number of pegRNAs generated in each codon of target gene.
    Also visualizes the identitity of mutant amino acid for SNPs in the input dataframe.
    
    FOR PROPER FUNCTION, REQUIRES USER INCLUSION OF HGVSp INFORMATION IN RF. See documentation for more info.
    
    Returns matrix_snp, matrix_ins, matrix_del, matrix_neutral, fig. First four return objects correspond with
    matrices used for plotting. Fig = figure, which can be saved/exported.
    
    Parameters
    ----------
    rf
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs generated by run().
        
    gene_name
        *type = str*
        
        Gene's Hugo_Symbol (i.e. name).
        
    start_end_cds
        *type = list*
        
        A 2-d list containing the start/end locations of each region of the coding sequence (CDS) for the gene's selected transcript.
        See documentation for example and precise specifications of format. 
        
    strand
        *type = str*
        
        Strand that transcript is on. Options are '+' or '-'.
        
    plot
        *type = boolean ('True' or 'False')*
        
        If true, show stacked bar plot. If false, simply return the matrices.
        

    """
    
    
    zero_array = lollipop_library(rf, gene_name, start_end_cds, strand, plot=False)
    
    df = rf
    
    amino_acids = ['ala','arg','asn','asp','cys','glu','gln','gly',
    'his','ile','leu','lys','met','phe','pro','ser',
    'thr','trp','tyr','val','TER']

    amino_acids = [i.upper() for i in amino_acids]
    dictionary = dict(zip(amino_acids, range(0,22)))
    
    
    
    #gene_idx = unique_genes.index(gene_name)
    #codons = codon_locs[gene_idx]
    

    codons = []
    for i in start_end_cds:
        for k in range(i[0], i[1]+1):
            codons.append(k)

    gene_codons = []
    if strand=='+':
        for i in range(0,len(codons),3):
            gene_codons.append([codons[i], codons[i+1], codons[i+2]])

    elif strand=='-':
        codons = codons[::-1]
        for i in range(0,len(codons),3):
            gene_codons.append([codons[i], codons[i+1], codons[i+2]])
    
    gene_codons = np.asarray(gene_codons)
    
    len_gene = len(gene_codons)+1
    
    matrix_snp = np.zeros((len(amino_acids), len_gene))
    matrix_neutral = np.zeros((1, len_gene))
    matrix_ins = np.zeros((1, len_gene))
    matrix_del = np.zeros((1, len_gene))
    
    
    
    
    for i, val in df.iterrows():
        
        mut_class = val['classification']
        
        typ = val['variant_type']
        
        #find codon
        i,j = np.where(gene_codons == val['start'])
        if len(i)>0:
            cod = i[0]
            
        else:
            cod='none'
        
        
        if typ=='SNP':
            
            HGVSp = val['HGVSp']
            if type(HGVSp)==str:
                #extracting relevant information from HGVSp
                if any(i.isdigit() for i in HGVSp)==True:
                    a = HGVSp.split('p.')[1]
                    split_a = re.split(r'(\d+)', a)
                    ref = split_a[0]
                    mut = split_a[2][0:3]
                    cod = int(split_a[1])

                 #if there's a valid HGVSp

                    mut_idx = dictionary[mut.upper()]
                    matrix_snp[mut_idx][int(cod)-1]+= 1#cod-1 corrects for shift by 1
                else:
                    continue
            else:
                continue #no HGVSp provided, continue
          
        
        if type(cod)!=str:
            if typ=='INS':
                matrix_ins[0][cod-1]+=1
            if typ=='DEL':
                matrix_del[0][cod-1]+=1

            if mut_class=='neutral':
                matrix_neutral[0][cod-1]+=1
        else:
            continue
                    
                
            
    #plotting
    if plot==True:
        fig, ax = plt.subplots(5,1,figsize=(25,20), sharex=True, gridspec_kw={ 'height_ratios':[3,8,0.7,.7,.7]})

        g = sns.heatmap(matrix_snp, cmap='Blues', ax=ax[1],xticklabels=10,alpha=1, cbar=False,cbar_kws={"orientation": "horizontal", "pad":0.05, "label":'Number of Occurences'})

        
        cmap_name='cyan1'
        colors = ['white','Grey']
        cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=4)
        h = sns.heatmap(matrix_neutral, cmap=cmap, ax=ax[2],xticklabels=10,alpha=1,cbar=False,cbar_kws={"orientation": "horizontal", "pad":0.05, "label":'Number of Occurences'})
        
        f = sns.heatmap(matrix_ins, cmap='Greens', vmin=0,ax=ax[3],xticklabels=10,alpha=1,cbar=False, cbar_kws={"orientation": "horizontal", "pad":0.05, "label":'Number of Occurences'})

        q = sns.heatmap(matrix_del, cmap='Purples', ax=ax[4],xticklabels=10,alpha=1,cbar=False,cbar_kws={"orientation": "horizontal", "pad":0.05, "label":'Number of Occurences'})
        
        zero = np.zeros(393)+9
        zero[0:41]=1 #AD1
        zero[42:91]=2 #AD2
        zero[101:291]=3 #DBD
        zero[304:321]=4 #NLS
        zero[325:355]=5 #TD
        zero[363:392]=6 #
        #p = sns.heatmap([zero], cmap='Pastel1', ax=ax[0], cbar=False, yticklabels = False,xticklabels = 20)#,cbar_kws = dict(use_gridspec=True,cax=ax[1][2]))
        
        #p.tick_params(bottom=False) 


        f.tick_params(bottom=False) 
        g.tick_params(bottom=False) 

        h.tick_params(bottom=False) 


                
        ax[0].bar(range(1, len_gene),zero_array[0], color='tab:blue',alpha=0.7, label='SNP')
        ax[0].bar(range(1, len_gene),zero_array[1], bottom=zero_array[0], color='tab:green',alpha=0.7, label='INS')
        ax[0].bar(range(1, len_gene),zero_array[2], bottom=zero_array[0]+zero_array[1], color='tab:purple', alpha=0.7,label='DEL')
        ax[0].bar(range(1, len_gene),zero_array[3],bottom=zero_array[0]+zero_array[1]+zero_array[2], color='tab:grey',alpha=0.7, label='ONP (neutral)')


        
        ax[0].set_xlim(1,len_gene)
        ax[0].set_ylabel('Number of pegRNAs', fontsize=17) 
        ax[4].set_xlabel('Codon', fontsize=15)
        ax[1].set_yticklabels(amino_acids, fontsize=15, rotation=0)
        ax[2].set_yticklabels(['NEUTRAL'], fontsize=15, rotation=0)
        ax[3].set_yticklabels(['INS'], fontsize=15, rotation=0)

        ax[4].set_yticklabels(['DEL'], fontsize=15, rotation=0)



        
        ax[0].set_title(gene_name, fontsize=25)

        for i in range(23):
            ax[1].axhline(i, color='white', lw=2)

        #h = sns.heatmap([zero_array], cmap='Reds', ax=ax[0])
        plt.subplots_adjust(wspace=0, hspace=0.05)
    
    #fig.colorbar(g, orientation="horizontal")


    
        return matrix_snp,matrix_ins, matrix_del, matrix_neutral, fig
    
    return matrix_snp,matrix_ins, matrix_del, matrix_neutral

#------------ClinVar variation id reader-----------------

def clinvar_VCF_translator(filepath, variation_ids):
    """
    Function that takes a clinvar.vcf.gz file containing information about variants, as well as a list of variation ID numbers,
    and returns a pandas dataframe containing the variants in a format that can be used by PEGG to design pegRNAs.
    
    Parameters
    ----------
    filepath
        *type = str*
        
        Filepath to the clinvar.vcf.gz file.
        
    variation_ids
        *type = list*
        
        List of variation IDs that the user wants to convert to .
        
    """
    gene = []
    ref = []
    alt = []
    chrom = []
    start = []
    end = []
    CLNHGVS = []
    var_type = []
    CLNDN = []
    CLNSIG = []
    allele_id = []
    
    var_id = []

    #for translating between Clinvar variation classes, and PEGG-readable versions...
    a = ['Deletion', 'Duplication', 'Indel', 'Insertion', 'Inversion',
        'Microsatellite', 'Variation', 'single_nucleotide_variant']
    b = ['DEL', 'Duplication', 'Indel', 'INS', 'Inversion',
        'Microsatellite', 'Variation', 'SNP']

    var_dict = dict(zip(a,b))


    for variant in VCF(filepath): # or VCF('some.bcf')
        if int(variant.ID) in variation_ids:
            
            var_id.append(int(variant.ID))
            allele_id.append(variant.INFO.get('ALLELEID'))
            #change it to make it flexible for a list
            gene.append(variant.INFO.get('GENEINFO'))
            CLNSIG.append(variant.INFO.get('CLNSIG'))
            var_type.append(variant.INFO.get('CLNVC'))
            CLNHGVS.append(variant.INFO.get('CLNHGVS'))
            CLNDN.append(variant.INFO.get('CLNDN'))
            ref.append(variant.REF)
            chrom.append(variant.CHROM)
            start.append(variant.start)
            end.append(variant.end)

            if len(variant.ALT)==1:
                alt.append(variant.ALT[0])
            elif len(variant.ALT)>1:
                tot = len(variant.ALT)-1
                alt.append(variant.ALT[0])
                for i in range(tot):
                    var_id.append(int(variant.ID))
                    allele_id.append(variant.INFO.get('ALLELEID'))
                    #change it to make it flexible for a list
                    gene.append(variant.INFO.get('GENEINFO'))
                    CLNSIG.append(variant.INFO.get('CLNSIG'))
                    var_type.append(variant.INFO.get('CLNVC'))
                    CLNHGVS.append(variant.INFO.get('CLNHGVS'))
                    CLNDN.append(variant.INFO.get('CLNDN'))
                    ref.append(variant.REF)
                    chrom.append(variant.CHROM)
                    start.append(variant.start)
                    end.append(variant.end)
                    alt.append(variant.ALT[i+1])
                    #d1.append(variant.ID)

    var_type = [var_dict[i] for i in var_type]
    gene = [i.split(':')[0] for i in gene] #extracting just gene information...

    #fixing issue where Clinvar variants start one position early for some reason...
    start = np.asarray(start)+1
    
    d1 = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 
          'Tumor_Seq_Allele2', 'Variant_Type', 'Variation_ID','Allele_ID','CLNSIG', 'CLNHGVS', 'CLNDN']

    combined = [gene, chrom, start, end, ref, alt, var_type, var_id, allele_id,
               CLNSIG, CLNHGVS, CLNDN]
    
    d = dict(zip(d1, combined))
    clinvar = pd.DataFrame(data = d)

    return clinvar


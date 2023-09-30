import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import Bio.Seq
import importlib

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


def neutral_substitutions(gene_name, chrom, strand, start_end_cds, chrom_dict):
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
        
    chrom_dict
        *type = dict*
        
        Dictionary containing the reference genome. See genome_loader() 
    
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
    
    #loading in reference genome from peg engine module
    seq1 = chrom_dict[chrom]

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

    g_cod = g_cod[g_cod['classification']=='neutral'].reset_index().drop(columns='index')

    #------and finally fixing the issue where SNPs are labelled as ONPS...

    for i, val in g_cod.iterrows():
        s = val['Start_Position']
        e = val['End_Position']

        r = val['Reference_Allele']
        a = val['Tumor_Seq_Allele2']

        diff_loc = []
        for index, val2 in enumerate(r):
            if val2!=a[index]:
                diff_loc.append(index)

        if len(diff_loc)==1:
            g_cod.loc[i, 'Variant_Type'] = 'SNP'
            g_cod.loc[i, 'Start_Position'] = s + diff_loc[0]
            g_cod.loc[i, 'End_Position'] = s + diff_loc[0]
            g_cod.loc[i, 'Reference_Allele'] = r[diff_loc[0]]
            g_cod.loc[i, 'Tumor_Seq_Allele2'] = a[diff_loc[0]]

    return g_cod


def safe_muts(num_muts, chrom_dict, organism='human'):
    """ 
    Function for generating a list of safe-targetting pegRNAs/gRNAs.
    NOTE: THIS ONLY WORKS FOR GRCh37 and GRCm38!!!
    Alternatively, generate these mutations/guides using these reference builds, and then generate your own guides
    targeting mutations in your desired reference build.

    This is a random subset of 100,000 "safe regions" taken from Morgens et al., 2017 (https://doi.org/10.1038/ncomms15178)

    Parameters
    -----------
    num_muts
        *type = int*
        
        Number of safe-targetting mutations to select/generate.

    chrom_dict
        *type = dict*

        Dictionary containing chromosomes in dictionary format for parsing by PEGG.

    organism
        *type = str*

        Choices = "human" or "mouse". Generates safe targeting guides in human or mouse genome.
    """

    assert organism in ['human', 'mouse'], "Pick a valid reference organism ('human' or 'mouse')"

    if organism == 'human':
        h = importlib.resources.files(__package__).joinpath('human_safe_regions_bassik.csv')
        d = pd.read_csv(h)

    if organism == 'mouse':
        m = importlib.resources.files(__package__).joinpath('mouse_safe_regions_bassik.csv')
        d = pd.read_csv(m)

    #dealing with mix of strings and ints in column because of pandas saving/loading issue with dtypes
    l = list(range(1,23))
    l = [str(i) for i in l]

    #choose regions without replacement
    a = list(range(len(d)))
    rand_choices = np.random.choice(a, size=num_muts, replace=False)
    d = d.iloc[rand_choices].reset_index()

    #now generate fake mutations (same REF and ALT) in cBioPortal format
    loc = []
    chrom_list = []
    dna_seq = []
    for i, val in d.iterrows():
        chrom = val['Chromosome']

        if chrom in l: #fixing int/str issue
            chrom = int(chrom)

        s = val['Start']
        e = val['End']
        #choose a random spot in the safe region as well...
        choice = np.random.randint(s,e)
        dna = chrom_dict[chrom][choice-1].upper() #-1 for 0 based indexing

        loc.append(choice)
        dna_seq.append(dna)
        chrom_list.append(chrom)


    col_labels = [ 'Chromosome', 'Variant_Type', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'classification']
    cols = [chrom_list, 'SNP', loc, loc, dna_seq, dna_seq, 'safe-targeting control']
    safe = pd.DataFrame(dict(zip(col_labels, cols)))

    safe = safe.sort_values(by=['Chromosome', 'Start_Position']).reset_index(drop=True)

    return safe


def aavs1_muts(chrom_dict, num_muts, genome_version = 'GRCh37'):
    """ 
    Function for generating list of "mutations" within AAVS1 locus as negative controls.
    Generates a list of SNPs that contain idential reference and alternate alleles.
    Intron 1-2 of PPP1R12C (AAVS1 locus): GRCh37 transcript = ENST00000263433.3 | GRCh38 transcript = ENST00000263433.8
    
    Not currently used in the pipeline, but provided anyway.

    Parameters
    -----------
    chrom_dict
        *type = dict*

        Dictionary containing chromosomes in dictionary format for parsing by PEGG.

    num_muts
        *type = int*
        
        Number of safe-targetting mutations to select/generate.

    genome_version
        *type = str*

        Options = 'GRCh37', 'GRCh38'

    """

    assert genome_version in ['GRCh37', 'GRCh38'], "Pick a valid reference genome ('GRCh37' or 'GRCh38')"

    if genome_version == 'GRCh37':
        #give 500 bp of buffer on either side of the intron
        start_intron = 55624164 + 500
        end_intron = 55628590 - 500

    elif genome_version == 'GRCh38':
        #give 500 bp of buffer on either side of the intron
        start_intron = 55112796 + 500
        end_intron = 55117222 - 500

    chr19 = chrom_dict[19]

    #select a random set of sequence locations in AAVS1
    snp_locs = np.random.random_integers(low = start_intron, high=end_intron, size=num_muts)

    #go through and generate the fake SNPs (same WT and REF allele) for these locations
    dna_seq = []
    for i in snp_locs:
        dna = chr19[int(i)-1]
        dna_seq.append(dna.upper())

    col_labels = ['Hugo_Symbol', 'Chromosome', 'Variant_Type', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'classification']
    cols = ['AAVS1', 19, 'SNP', snp_locs, snp_locs, dna_seq, dna_seq, 'AAVS1-targeting control']
    aavs1 = pd.DataFrame(dict(zip(col_labels, cols)))

    return aavs1


def nontargeting_guides(num_guides, edit_type='prime'):
    """ 
    Function for generating a list of non-targetting guides (in the human genome).
    Future versions will expand to mouse as well...Unprocessed files for mouse located within this module.
    Guides taken from: https://www.nature.com/articles/nmeth.4423 (https://doi.org/10.1038/nmeth.4423)

    Parameters
    -----------
    num_guides
        *type = int*

        Number of guides to non-targetting generate. Max = 1000.

    edit_type
        *type = str*

        Options = 'base', 'prime'
    """

    assert edit_type in ['base', 'prime'], "Choose valid edit_type (choices = 'prime', 'base')"

    if num_guides>1000:
        print('Only 1000 non-targeting guides available; will add 1000 guides.')
        num_guides = 1000

    if edit_type == 'prime':
        pr = importlib.resources.files(__package__).joinpath('non_targeting_human_prime.csv')
        df = pd.read_csv(pr)

    if edit_type == 'base':
        ba = importlib.resources.files(__package__).joinpath('non_targeting_human_base.csv')
        df = pd.read_csv(ba)

    df = df[:num_guides]

    return df
    


def library_maker(mutant_input, gene_name, chrom_dict, fraction_safetarget=0.05,organism='human', 
                  fraction_silent=0, chrom=None, strand=None, start_end_cds=None):
    """ 
    Compiles the different library design functions to aggregate all of the variants for a particular gene,
    include silent substitution controls at a designated %, and non-targeting controls. This can then be fed
    into the pegg2.run() or base_editing.run_base() functions.

    NOTE: THIS ONLY WORKS FOR GRCh37 and GRCm38!!!
    Alternatively, generate these mutations/guides using these reference builds, and then generate your own guides
    targeting mutations in your desired reference build.

    Parameters
    -----------
    mutant_input
        *type = pd.DataFrame*

        DataFrame containing all of the input mutations. Doesn't require filtration for gene of interest.
    
    gene_name
        *type = str*

        Gene name to select from the mutant_input dataframe. 

    chrom_dict
        *type = dict*

        Dictionary containing the reference genome. See genome_loader() 

    fraction_safetarget
        *type = float*

        Value from 0 to 1 that corresponds to the fraction of safe-targetting mutations to include in the library.
        (i.e. what fraction of the mutant_input)

    organism
        *type = str*

        Options = 'human' or 'mouse'. Determines which of the organisms to generate safe targeting guides for.

    fraction_silent
        *type = float*

        Value from 0 to 1 that corresponds to the fraction of silent mutations to include in the library.
        (i.e. what fraction of the mutant_input)

    chrom
        *type = int or str*

        Chromosome that the gene of interest falls on.

    strand
        *type = str*

        '+' or '-' -- corresponds with which strand the transcript falls on.

    start_end_cds
        *type = list*

        Nested list that contains the CDS locations of the gene transcript of interest, in the + strand orientation.
    
    """
    assert (fraction_silent>=0) & (fraction_silent<=1), "Choose valid fraction_silent (between 0 and 1)"
    assert (fraction_safetarget>=0) & (fraction_safetarget<=1), "Choose valid fraction_nontarget (between 0 and 1)"
    assert organism in ['human', 'mouse'], "Pick a valid reference organism ('human' or 'mouse')"

    #aggregate variants corresponding to gene
    agg_muts = mutation_aggregator(mutant_input, gene_name)
    agg_muts['classification'] = 'variant'
    total_size = len(agg_muts)

    #generate neutral/silent substitution controls
    if fraction_silent>0:
        neutrals = neutral_substitutions(gene_name, chrom, strand, start_end_cds, chrom_dict)

        num_silent = min(int(total_size*fraction_silent), len(neutrals))

        a = list(range(len(neutrals)))
        rand_choices = np.random.choice(a, size=num_silent, replace=False)
        neutrals = neutrals.iloc[rand_choices].reset_index()

        #generate safe-targetting variants
        num_safe = int(total_size*fraction_safetarget)
        safe = safe_muts(num_safe, chrom_dict, organism)

        return pd.concat([agg_muts, neutrals, safe]).reset_index(drop=True)

    else:
        #generate safe-targetting variants
        num_safe = int(total_size*fraction_safetarget)
        safe = safe_muts(num_safe, chrom_dict, organism)

        return pd.concat([agg_muts, safe]).reset_index(drop=True)

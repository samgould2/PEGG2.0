import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
plt.rc('font', family='Helvetica')
import Bio.Seq
import gzip
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.pairwise2 import format_alignment
import warnings
import regex as re
import pickle
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from cyvcf2 import VCF
import crisporEffScores #python script included in directory; from CRISPOR website github with modifications

warnings.filterwarnings('ignore')

#--- GLOBAL VARIABLES ---------
MATCH_SCORE = 1
MISMATCH_SCORE = -0.5
OPEN_GAP_SCORE = -5
EXTEND_GAP_SCORE = -0.1
TARGET_END_GAP_SCORE = 0
QUERY_END_GAP_SCORE = 0

def make_aligner():
    """ 
    Aligner function from Bio.Align.PairwiseAligner with custom parameters.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = MATCH_SCORE
    aligner.mismatch_score = MISMATCH_SCORE
    aligner.open_gap_score = OPEN_GAP_SCORE
    aligner.extend_gap_score = EXTEND_GAP_SCORE
    aligner.target_end_gap_score = TARGET_END_GAP_SCORE
    aligner.query_end_gap_score = QUERY_END_GAP_SCORE
    return aligner

aligner = make_aligner()

class mutation:
    """ 
    Class for storing information about individual mutations. Used by functions throughout.
    Prevents the use of disorganized lists.
    """
    def __init__(self, wt_w_context, alt_w_context, left_seq, right_seq, var_type, ref_seq, alt_seq, chrom=None, genome=None):
        """ 
        change seq_start and seq_end to correspond with chromosome coordinates
        or some other organizational method...
        """
        
        #WT sequence
        self.wt_forward = wt_w_context.upper()
        rc = str(Bio.Seq.Seq(self.wt_forward).reverse_complement())
        self.wt_rc = rc

        #ALT sequence
        self.mut_forward = alt_w_context.upper()
        rc_mut = str(Bio.Seq.Seq(self.mut_forward).reverse_complement())
        self.mut_rc = rc_mut

        #seq start refers to the WIDE sequence start (i.e. the WT sequence with e.g. 60 nt of context)
        self.seq_start = 0
        self.seq_end = len(self.wt_forward)
        self.chrom = chrom
        self.genome_build = genome

        #and splitting it up into the constituent components
        self.left_seq = left_seq.upper()
        self.right_seq = right_seq.upper()
        self.ref_seq = ref_seq.upper()
        self.alt_seq = alt_seq.upper()

        #and do the same to get it into the RC format

        #RIGHT/LEFT ARE SWAPPED ON REVERSE COMPLEMENT
        self.left_seq_rc = str(Bio.Seq.Seq(right_seq).reverse_complement())
        self.right_seq_rc = str(Bio.Seq.Seq(left_seq).reverse_complement())
        self.ref_seq_rc = str(Bio.Seq.Seq(ref_seq).reverse_complement())
        self.alt_seq_rc = str(Bio.Seq.Seq(alt_seq).reverse_complement())

        self.variant_type = var_type
        #self.alt_size = (positive or negative integer)
        
        #information about the PAM sequence location
        self.PAM_idx_forward = None
        self.PAM_idx_rc = None


#-----------functions------------
def genome_loader(filepath_gz):
    """
    Takes in filepath of human or mouse genome and returns dictionary of chromosome sequences that PEGG can parse.
    Tested only on human and mouse genomes.

    Returns (1) chromosome dictionary, and (2) file names of the chromosomes that are stored (for manual checking of errors).
    
    Parameters
    -----------
    filepath_gz
        *type = str*
        
        The filepath to the .gz file holding the reference genome file.
 
    """
    #------loading in reference genome and organizing it into a 2-d list by chromosome---------------------
    wrong = ["alternate", "unplaced", "unlocalized", "patch", "mitochondrion"]
    filtered = []
    chrom_list = []
    seq_list = []

    with gzip.open(filepath_gz, "rt") as handle:
        #records = list(SeqIO.parse(handle, "fasta")) #about 4 Gb in  memory

        for i in SeqIO.parse(handle, "fasta"):
            ii = i.description
            ignore=False 

            for key in wrong:
                if key in ii:
                    ignore = True
            
            if ignore==False:

                x = ii.find('chromosome')
                chrom = ii[x:].split(' ')[1]
                chrom = chrom.split(',')[0]
                if chrom not in ['X', 'Y']:
                    chrom = int(chrom)
                chrom_list.append(chrom)
                seq_list.append(i.seq)

                filtered.append(ii)

  
    return dict(zip(chrom_list, seq_list)), filtered

def df_formatter(df, chrom_dict, context_size = 120): 
                 #cols_to_save = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Consequence', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2',  'HGVSc', 'HGVSp_Short']):

    """ 
    Takes in variants (in cBioPortal format!) and outputs dataframe with REF and ALT oligos with designated context_size
    that can be used by PEGG for creating pegRNAs.

    Parameters
    -----------
    df
        *type = pd.DataFrame*
        
        The dataframe of input variants in cBioPortal format.

    chrom_dict
        *type = dict*

        Dictionary generated by genome_loader() that holds the chromosome sequences.

    context_size
        *type = int*

        The amount of context/flanking sequence on either side of the mutation to generate. For larger variants, set this larger.
        Default = 120. e.g. AAA(A/G)AAA = context_size of 3
 
    """

    wt_w_context = []
    alt_w_context = []
    left_context_list = []
    right_context_list = []
    ref_allele = []
    alt_allele = []

    seq_start = []
    seq_end = []

    to_drop = []

    for i, val in df.iterrows():
        vt = val['Variant_Type']
        s = val['Start_Position']
        e = val['End_Position']
        ref = val['Reference_Allele']
        alt = val['Tumor_Seq_Allele2']
        chrom = val['Chromosome']

        if chrom not in ['X','Y']:
            chrom = int(chrom)
       
        chr_seq = chrom_dict[chrom].upper()

        if vt in ['SNP', 'ONP', 'DNP', 'INDEL']:
            ref = ref
            alt = alt
            #assert ref == chr_seq[s-1:e], print(ref, chr_seq[s-1:e])
            left_context = chr_seq[s-1-context_size:s-1]
            right_context = chr_seq[e:e+context_size]

        elif vt =='INS':
            ref = ''
            alt = alt
            #left_context = chr_seq[s-1-context_size:s+1] #need to do this since INS reference alleles are blank
            left_context = chr_seq[s-1-context_size:s]
            right_context = chr_seq[e-1:e+context_size]

        elif vt=='DEL':
            ref = ref
            alt = ''
            left_context = chr_seq[s-1-context_size:s-1]
            right_context = chr_seq[e:e+context_size]

        wt_seq = left_context + ref + right_context
        alt_seq = left_context + alt + right_context

        

        start = s-context_size
        end = e+context_size

        if str(chr_seq[start-1:end])!=str(wt_seq): #changed from assert to allow it to continue to run

            print(f"Error in mutant # {i}; Dropped from pegRNA generation\nVariant Type = {vt}| REF = {ref} | ALT = {alt}\n{chr_seq[start-1:end]}\n{str(wt_seq)}")
            to_drop.append(i)
        else:
            seq_start.append(start)
            seq_end.append(end)
            
            wt_w_context.append(str(wt_seq))
            alt_w_context.append(str(alt_seq))
            left_context_list.append(str(left_context))
            right_context_list.append(str(right_context))
            ref_allele.append(str(ref))
            alt_allele.append(str(alt))

    cols_to_save = list(df.keys())

    df_new = df[cols_to_save]
    df_new = df_new.drop(index=to_drop).reset_index(drop=True)

    df_new['seq_start'] = seq_start
    df_new['seq_end'] = seq_end
    df_new['wt_w_context'] = wt_w_context
    df_new['alt_w_context'] = alt_w_context
    df_new['left_context'] = left_context_list
    df_new['right_context'] = right_context_list
    df_new['REF'] = ref_allele
    df_new['ALT'] = alt_allele
    df_new = df_new.reset_index()
    
    return df_new

def mut_formatter(wt, alt):
    """ 
    Formats mutations for pegRNA generation. Takes in WT and ALT sequences, returns necessary parameters for gRNA/pegRNA generation.
    Note: This function will break with INDELs and complex variants since it relies on a simple alignment.

    Parameters
    -----------
    wt
        *type = str*
        
        WT sequence

    alt
        *type = str*

        Mutant/alternate sequence.
    """

    alignments = aligner.align(wt, alt)

    WT_align = alignments[0].aligned[0]
    alt_align = alignments[0].aligned[1]

    if len(WT_align) == 1:
        #Substitution (SNP or ONP)

        #find indexes of non-matching values to get right/left context sequences and mutant sequence
        non_matching_idx = []

        for i, val in enumerate(list(zip(wt, alt))):
            if val[0] != val[1]:
                non_matching_idx.append(i)
        
        assert len(non_matching_idx) !=0, f'No mutation detected in alt sequence provided | row index = {i}'

        if len(non_matching_idx)==1:
            var_type = 'SNP'
        else:
            var_type = 'ONP'

        left_seq = wt[:non_matching_idx[0]]
        ref_seq = wt[non_matching_idx[0]:non_matching_idx[-1]+1]
        alt_seq = alt[non_matching_idx[0]:non_matching_idx[-1]+1]
        right_seq = wt[non_matching_idx[-1]+1:]


    elif len(WT_align)>1:
        if len(alt)>len(wt):
            var_type = 'INS'
        elif len(alt)<len(wt):
            var_type='DEL'

        left_seq = wt[WT_align[0][0]:WT_align[0][1]]
        right_seq = wt[WT_align[1][0]:WT_align[1][1]]

        if var_type=='INS':
        #INS or DEL
            ref_seq = ''
            alt_seq = alt[alt_align[0][1]:alt_align[1][0]]
            
        elif var_type =='DEL':
            ref_seq = wt[WT_align[0][1]:WT_align[1][0]]
            alt_seq = ''
            

    return var_type, left_seq, ref_seq, alt_seq, right_seq


def primedesign_formatter(seq):
     
    """ 
    Takes as input sequence in prime design format e.g. AATTCCG(G/C)AATTCGCT. Returns necessary parameters for pegRNA generation.

    Parameters
    -----------
    seq
        *type = str*
        
        The PrimeDesign formatted sequence.
    """

    start = seq.find("(")

    end = seq.find(")")

    replace_seq = seq[start:end+1]
    #if '/' in replace_seq:
        #if '+' in replace_seq:
            #throw an error

    loc_replace = replace_seq.find('/')
    ref_seq = replace_seq[1:loc_replace]
    alt_seq = replace_seq[loc_replace+1:-1]

    left_seq = seq[:start]
    right_seq = seq[end+1:]

    wt_w_context = left_seq + ref_seq + right_seq
    alt_w_context = left_seq + alt_seq + right_seq

    if len(ref_seq)==0:
        var_type='INS'

    else:
        if len(ref_seq)==len(alt_seq):
            if len(ref_seq)==1:
                var_type='SNP'
            else:
                var_type='ONP'

        else:
            if len(alt_seq) == 0:
                var_type = 'DEL'
            else:
                var_type = 'INDEL'

    return var_type, wt_w_context, alt_w_context, left_seq, right_seq, ref_seq, alt_seq


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
    b = ['DEL', 'Duplication', 'INDEL', 'INS', 'Inversion',
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

def PAM_finder(seq, PAM):
    """ 
    Finds indeces of PAM sequences in sequence.

    Parameters
    ------------
    seq
        *type = str*

        Sequence to search for PAM sequences.

    PAM
        *type = str*
        
        PAM sequence for searching. Formatting: e.g. "NGG"

        N = [A|T|C|G]
        R = [A|G]
        Y = [C|T]
        S = [G|C]
        W = [A|T]
        K = [G|T]
        M = [A|C]
        B = [C|G|T]
        D = [A|G|T]
        H = [A|C|T]
        V = [A|C|G]
    """
    PAM_dict = {"A":"A",
                "G":"G",
                "C":"C",
                "T":"T",
                "N":"[A|T|C|G]",
                "R":"[A|G]",
                "Y":"[C|T]",
                "S":"[G|C]",
                "W":"[A|T]",
                "K":"[G|T]",
                "M":"[A|C]",
                "B":"[C|G|T]",
                "D":"[A|G|T]",
                "H":"[A|C|T]",
                "V":"[A|C|G]"}

    PAM_new = ""
    for i in PAM:
        PAM_new += PAM_dict[i]

    p = re.compile(PAM_new)

    iterator = p.finditer(seq, overlapped=True)
    PAM_holder = []
    for match in iterator:
        PAM_holder.append(match.span())

    return PAM_holder


def eligible_PAM_finder(mut, PAM, max_RTT_length, proto_size=19):
    """ 
    Determines eligible PAM sequences for creating pegRNAs.
    Returns PAM sequence locations on (1) forward, and (2) reverse-complement strand

    Parameters
    ------------
    mut
        *type = mutation class*

        See class: mutation

    PAM
        *type = str*
        
        PAM sequence for searching. Formatting: e.g. "NGG"

        N = [A|T|C|G]
        R = [A|G]
        Y = [C|T]
        S = [G|C]
        W = [A|T]
        K = [G|T]
        M = [A|C]
        B = [C|G|T]
        D = [A|G|T]
        H = [A|C|T]
        V = [A|C|G]

    max_RTT_length
        *type = int*

        Max RTT length for pegRNAs being designed.

    proto_size
        *type = int*

        Size of protospacer being used. Default = 19 (G+19).
    """

    forward = PAM_finder(mut.wt_forward, PAM)
    reverse = PAM_finder(mut.wt_rc, PAM)

    var_type = mut.variant_type

    #first do the forward sequence
    left_len_F = len(mut.left_seq)
    left_len_R = len(mut.left_seq_rc)

    if var_type =='INS':
        alt_size = len(mut.alt_seq)
        m_start_F = left_len_F - (max_RTT_length - alt_size - 3)
        m_end_F = left_len_F + 4
        
        m_start_R = left_len_R - (max_RTT_length - alt_size - 3)
        m_end_R = left_len_R + 4

    else:
        alt_size = len(mut.alt_seq)
        m_start_F = left_len_F - (max_RTT_length - alt_size - 3)
        m_end_F = left_len_F + 3

        m_start_R = left_len_R - (max_RTT_length - alt_size - 3)
        m_end_R = left_len_R + 3


    eligible_PAM_start_F = max(m_start_F, proto_size + 3) #
    eligible_PAM_end_F = m_end_F
    
    eligible_PAM_start_R = max(m_start_R, proto_size + 3) #
    eligible_PAM_end_R = m_end_R

    eligible_PAMS_F = []
    eligible_PAMS_R = []

    for i in forward:
        #check start of PAM site
        s = i[0]
        if (s>=eligible_PAM_start_F) and (s<=eligible_PAM_end_F):
            eligible_PAMS_F.append(i)

    for i in reverse:
        #check start of PAM site
        s = i[0]
        if (s>=eligible_PAM_start_R) and (s<=eligible_PAM_end_R):
            eligible_PAMS_R.append(i)

    return eligible_PAMS_F, eligible_PAMS_R

def pegRNA_generator(mut, PAM, orientation, proto_size, RTT_lengths, PBS_lengths):
    """ 
    Generates pegRNAs for a given mutation.

    Parameters
    ------------
    mut
        *type = mutation class*

        See class: mutation

    PAM
        *type = str*
        
        PAM sequence for searching. Formatting: e.g. "NGG"

        N = [A|T|C|G]
        R = [A|G]
        Y = [C|T]
        S = [G|C]
        W = [A|T]
        K = [G|T]
        M = [A|C]
        B = [C|G|T]
        D = [A|G|T]
        H = [A|C|T]
        V = [A|C|G]

    orientation
        *type = str*

        '+' or '-'

    proto_size
        *type = int*

        Size of protospacer being used. Default = 19 (G+19).

    RTT_lengths
        *type = list*

        List containing RTT lengths to design pegRNAs for.

    PBS_length
        *type = list*

        List containing PBS lengths to desing pegRNAs for.
    """


    #------function----
    protospacer_seqs = []
    protospacer_wide_seqs = []
    PAM_start_list = []
    PAM_end_list = []
    PAM_list = []

    RTT_seqs = []
    RTT_size = []
    PBS_seqs = []
    PBS_size = []
    RTT_PBS_seqs = []

    #descriptors of pegRNAs...
    distance_to_nick = [] #distance from the start of the mutation to the nick introduced by nCas9 (-3 from PAM)
    RHA = [] #right homology arm (distance from end of mutation to end of RTT)
    PAM_disrupted_list = []
    proto_disrupted_list = []

    #Get mutation info, depending on whether we're looking on the + or the - strand
    if orientation == '+':
        left_len = len(mut.left_seq)
        ref_len = len(mut.ref_seq)
        alt_len = len(mut.alt_seq)
        alt_seq = mut.alt_seq
        seq_F = mut.wt_forward

        PAM_seqs = mut.PAM_idx_forward

    elif orientation == '-':
        left_len = len(mut.left_seq_rc)
        ref_len = len(mut.ref_seq_rc)
        alt_len = len(mut.alt_seq_rc)
        alt_seq = mut.alt_seq_rc
        seq_F = mut.wt_rc

        PAM_seqs = mut.PAM_idx_rc

    #iterate over the PAM sequences:
    for i in PAM_seqs:
        PAM_start = i[0]
        PAM_end = i[1]

        RTT_start = PAM_start - 3

        left_RTT = seq_F[RTT_start:left_len]

        for RTT_length in RTT_lengths:

            remaining_length = RTT_length - (len(left_RTT) + alt_len)

            #make sure RTT is long enough
            #and make sure that there's sufficient sequence to actually pull from...
            if (remaining_length <0) or ((left_len+ref_len + remaining_length) > len(seq_F)):
                RTT_seqs.append(None)
                distance_to_nick.append(None)
                RHA.append(None)
                PAM_disrupted_list.append(None)
                proto_disrupted_list.append(None)
                PAM_start_list.append(None)
                PAM_end_list.append(None)
                PAM_list.append(None)
                protospacer_seqs.append(None)
                protospacer_wide_seqs.append(None)
                PBS_seqs.append(None)
                RTT_PBS_seqs.append(None)
                RTT_size.append(None)
                PBS_size.append(None)

            else: #able to design the pegRNA

                right_RTT = seq_F[left_len+ref_len:left_len+ref_len + remaining_length]
                RTT = left_RTT + alt_seq + right_RTT

                #and then determine if the PAM sequence is disrupted
                #protospacer disrupted would just be if the distance to nick is < 3
                #fish out the PAM sequence from RTT
                PAM_new = RTT[3:3+len(PAM)]
                pam_list = PAM_finder(PAM_new, PAM)

                PAM_disrupted = False
                if len(pam_list)==0:
                    PAM_disrupted = True

                proto_disrupted = False
                if len(left_RTT) < 3:
                    proto_disrupted = True

                #and finally reverse complement the RTT
                RTT = str(Bio.Seq.Seq(RTT).reverse_complement())

                #pull out the PAM sequence and protospacer as well
                PAM_sequence = seq_F[PAM_start:PAM_end]
                protospacer =  'G' + seq_F[PAM_start - proto_size:PAM_start]

                #extract wide protospacer for Rule set 2 prediction...
                protospacer_wide = seq_F[PAM_start - 20 - 4 : PAM_start + 3 + 3]
                
                #iterate over PBS lengths for adding things to list
                for PBS_length in PBS_lengths:
                    RTT_seqs.append(RTT)
                    RTT_size.append(RTT_length)
                    PBS_size.append(PBS_length)
                    distance_to_nick.append(len(left_RTT))
                    RHA.append(len(right_RTT))
                    PAM_disrupted_list.append(PAM_disrupted)
                    proto_disrupted_list.append(proto_disrupted)

                    #------add in PAM info, protospacer, + PBS/RTT_PBS
                    #record info about PAM sequences
                    if orientation == '+':
                        PAM_start_list.append(PAM_start)
                        PAM_end_list.append(PAM_end)
                    elif orientation == '-':
                        PAM_start_list.append(len(seq_F)-PAM_start)
                        PAM_end_list.append(len(seq_F)-PAM_end)

                    PAM_list.append(PAM_sequence)
                    protospacer_seqs.append(protospacer)
                    protospacer_wide_seqs.append(protospacer_wide)

                    #----and finally design the PBS sequences----
                    #NEED TO HAVE AN ERROR MESSAGE FOR PBS SEQUENCES LARGER THAN 17
                    PBS = seq_F[RTT_start-PBS_length:RTT_start]
                    PBS = str(Bio.Seq.Seq(PBS).reverse_complement())
                    PBS_seqs.append(PBS)

                    #---and then add it to the RTT to make the full 3' extension
                    RTT_PBS_seqs.append(RTT+PBS)

    cols = [PAM_start_list, PAM_end_list, PAM_list, orientation, protospacer_wide_seqs, protospacer_seqs, RTT_seqs, RTT_size, PBS_seqs, PBS_size, RTT_PBS_seqs, distance_to_nick, RHA, PAM_disrupted_list, proto_disrupted_list]
    col_labels = ["PAM_start", "PAM_end", "PAM", "PAM_strand","Protospacer_30", "Protospacer", "RTT", "RTT_length", "PBS", "PBS_length", "RTT_PBS", "Distance_to_nick", "RHA_size", "PAM_disrupted", "Proto_disrupted"]
    dtypes = ['int', 'int', 'str', 'str', 'str', 'str', 'str', 'int', 'str', 'int', 'str', 'int', 'int', 'bool', 'bool']

    dtype_dict = dict(zip(col_labels, dtypes))

    peg_df = pd.DataFrame(dict(zip(col_labels, cols))).dropna().reset_index().drop(columns='index').astype(dtype_dict)

    return peg_df

def sensor_generator(df, proto_size, before_proto_context=5, sensor_length=60, sensor_orientation = 'reverse-complement'):
    """ 
    Generates sensor sequence for quantification of pegRNA editing outcomes.
    Automatically puts sensor in reverse complement orientation with respect to protospacer.
    This is highly reccomended to reduce recombination during cloning/library preparation.

    Parameters
    -----------
    df
        *type = pd.DataFrame*

        Dataframe containing pegRNAs.

    proto_size
        *type = int*

        Size of protospacer. 

    before_proto_context
        *type = int*

        Default = 5. Amount of nucleotide context to put before the protospacer in the sensor

    sensor_length
        *type = int*

        Total length of the sensor in nt. Default = 60. 

    sensor_orientation
        *type = str*

        Options for sensor_orientation = 'reverse-complement' or'forward'.

    """

    assert sensor_orientation in ['reverse-complement', 'forward'], "Choose an eligible sensor orientation option ('reverse_comp' (default) or 'forward')"

    sensor_error = []
    sensor_wt_seqs = []
    sensor_alt_seqs = []

    for i, val in df.iterrows():
        wt = val['wt_w_context']
        alt = val['alt_w_context']
        ref_minus_alt = len(val['REF']) - len(val['ALT'])
        strand = val['PAM_strand']
        pam_start = val['PAM_start']
        pam_end = val['PAM_end']
        RTT_PBS = val['RTT_PBS']
        RHA = val['RHA_size']

        #if PAM is on negative strand; flip it
        if strand == '-':
            pam_start = len(wt) - pam_start
            pam_end = len(wt) - pam_end
            wt = str(Bio.Seq.Seq(wt).reverse_complement())
            alt = str(Bio.Seq.Seq(alt).reverse_complement())

        s1 = pam_start-(proto_size+1)-before_proto_context
        if s1 <0:
            sensor_wt_seqs.append(None)
            sensor_alt_seqs.append(None)
            sensor_error.append('insufficient context sequence provided; try increasing context size')

        else:
            proto_pam = wt[s1:pam_end]
            proto_pam_alt = alt[s1:pam_end]

            remain_len = sensor_length - len(proto_pam)
            if remain_len < 0:
                sensor_wt_seqs.append(None)
                sensor_alt_seqs.append(None)
                sensor_error.append('sensor too small; increase sensor size in parameters; or decrease before_proto_context')

            else:
                remain_len_alt = sensor_length - len(proto_pam) - ref_minus_alt

                e1 = pam_end + remain_len
                e2 = pam_end + remain_len_alt
                if e1>len(wt):
                    sensor_wt_seqs.append(None)
                    sensor_alt_seqs.append(None)
                    sensor_error.append('insufficient context sequence provided; try increasing context size')

                else:
                    rtt_region = wt[pam_end: e1]
                    rtt_region_alt = alt[pam_end:e2]

                    sensor_wt = proto_pam + rtt_region
                    sensor_alt = proto_pam_alt + rtt_region_alt

                    #AND FINALLY either keep it as is, or take the reverse complement
                    if sensor_orientation == 'reverse-complement':
                        sensor_wt = str(Bio.Seq.Seq(sensor_wt).reverse_complement())
                        sensor_alt = str(Bio.Seq.Seq(sensor_alt).reverse_complement())

                        rha_seq = RTT_PBS[:RHA]

                        #check that the homology overhang fits within the sensor sequence
                        if rha_seq not in sensor_wt:
                            sensor_wt_seqs.append(None)
                            sensor_alt_seqs.append(None)
                            sensor_error.append("sensor too small; increase sensor size in parameters; or decrease before_proto_context")

                        else:
                            sensor_wt_seqs.append(sensor_wt)
                            sensor_alt_seqs.append(sensor_alt)
                            sensor_error.append("No Error")
                    
                    if sensor_orientation == 'forward':
                        rha_seq = str(Bio.Seq.Seq(RTT_PBS[:RHA]).reverse_complement())

                        if rha_seq not in sensor_wt:
                            sensor_wt_seqs.append(None)
                            sensor_alt_seqs.append(None)
                            sensor_error.append("sensor too small; increase sensor size in parameters; or decrease before_proto_context")
                    
                        else:
                            sensor_wt_seqs.append(sensor_wt)
                            sensor_alt_seqs.append(sensor_alt)
                            sensor_error.append("No Error")


    df['sensor_wt'] = sensor_wt_seqs
    df['sensor_alt'] = sensor_alt_seqs
    df['sensor_orientation'] = sensor_orientation
    df['sensor_error'] = sensor_error

    return df

def ontarget_score(df):
    """ 
    Calls to crisporEffScores (taken from CRISPOR github) to generate on-target scores using Rule Set 2.
    """

    uniq_protos = list(np.unique(df['Protospacer_30']))

    #make sure they're all the right length
    uniq_protos_true = []
    wrong_length = []
    for val in uniq_protos:
        if len(val)!=30:
            wrong_length.append(val)
        else:
            uniq_protos_true.append(val)

    scores = crisporEffScores.calcAziScore(uniq_protos_true)

    dict_scores = dict(zip(uniq_protos_true, scores))
    dict_none = dict(zip(wrong_length, [None]*len(wrong_length)))

    dict_scores.update(dict_none)

    for protospacer in dict_scores.keys():
        df.loc[df["Protospacer_30"]== protospacer, "OnTarget_Azimuth_Score"] = dict_scores[protospacer]

    #for protospacer in dict_none.keys():
    #    df.loc[df["Protospacer_30"]== protospacer, "OnTarget_Azimuth_Score"] = None

    return df


def input_formatter(input_df, input_format, chrom_dict, context_size):
    """ 
    Master function for putting input mutations in the correct format for pegRNA/gRNA design.

    input_df
        *type = pd.DataFrame*

        Pandas dataframe that contains the input mutations.

    input_format
        *type = str*

        Options = 'cBioPortal', 'WT_ALT', 'PrimeDesign'. For 'WT_ALT', make sure to put headers as "WT" and "ALT". 
        For "PrimeDesign", put the header as "SEQ".

    chrom_dict
        *type = dict*

        Dictionary generated by genome_loader() that holds the chromosome sequences.

    context_size
        *type = int*

        The amount of context/flanking sequence on either side of the mutation to generate. For larger variants, set this larger.
        Default = 120. e.g. AAA(A/G)AAA = context_size of 3

    """

    assert input_format in ['cBioPortal', 'WT_ALT', 'PrimeDesign'], "input_format is not correct. Choose from ['cBioPortal', 'WT_ALT', 'PrimeDesign']."

    #format the input dataframe appropriately...
    if input_format == 'cBioPortal':
        #add assert statement to make sure all necessary info is included
        #MAKE COLS_TO_SAVE MORE FLEXIBLE??? to avoid weird errors...

        input_df = df_formatter(input_df, chrom_dict, context_size)
        input_df['mutation_idx'] = list(range(len(input_df)))

    elif input_format == 'WT_ALT':

        wt = []
        mut = []
        vts = []
        l = []
        ref = []
        alt =[]
        r = []
        for i, val in input_df.iterrows():

            var_type, left_seq, ref_seq, alt_seq, right_seq = mut_formatter(val['WT'], val['ALT'])
            vts.append(var_type)
            l.append(left_seq)
            ref.append(ref_seq)
            alt.append(alt_seq)
            r.append(right_seq)
            wt.append(val['WT'])
            mut.append(val['ALT'])

        mut_idx = list(range(len(wt)))
        cols = [mut_idx, wt, mut, vts, ref, alt, l, r]
        col_labels = ["mutation_idx", 'wt_w_context', 'alt_w_context', 'Variant_Type', 'REF', 'ALT', 'left_context', 'right_context']
        input_df = pd.DataFrame(dict(zip(col_labels, cols)))        

    elif input_format == 'PrimeDesign':

        wt = []
        mut = []
        vts = []
        l = []
        ref = []
        alt =[]
        r = []
        for i, val in input_df.iterrows():
            var_type, wt_w_context, alt_w_context, left_seq, right_seq, ref_seq, alt_seq = primedesign_formatter(val['SEQ'])
            vts.append(var_type)
            l.append(left_seq)
            ref.append(ref_seq)
            alt.append(alt_seq)
            r.append(right_seq)
            wt.append(wt_w_context)
            mut.append(alt_w_context)

        mut_idx = list(range(len(wt)))
        cols = [mut_idx, wt, mut, vts, ref, alt, l, r]
        col_labels = ["mutation_idx", 'wt_w_context', 'alt_w_context', 'Variant_Type', 'REF', 'ALT', 'left_context', 'right_context']
        input_df = pd.DataFrame(dict(zip(col_labels, cols)))     

    return input_df

def other_filtration(pegRNA_df, RE_sites=None, polyT_threshold=4):

    """
    Determines whether pegRNAs contain polyT termination sites or Restriction Enzyme Sites. Returns pegRNA_df with properties added.
    Note: RE_sites requires only A, T, C, or G bases (Ns will not work currently; need to do custom filtration)...
    These properties are checking for:
    
    (1) u6 terminator (polyT sequence)
        
    (2) RE site presence 
        
    Parameters
    ----------
    pegRNA_df
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()

    RE_sites
        *type = list or None*

        A list containing the RE recognition sites to filter (e.g. ['CGTCTC', 'GAATTC'] for Esp3I and EcoRI). Default = None (no filtration).

    polyT_threshold
        *type = int*

        The length of the polyT sequence to classify as a terminator. Default = 4.
    
    """
    #checking the type of dataframe

    prime=False
    base=False
    sensor=False

    df_keys = pegRNA_df.keys()
    if "RTT_PBS" in df_keys:
        prime=True
    else:
        base=True
    if 'sensor_wt' in df_keys:
        sensor=True

    terminator_sequence = 'T'*polyT_threshold

    contains_terminator = []
    RE_site = []

    if RE_sites != None:

        assert type(RE_sites)==list, "RE_sites must be in list format (e.g. ['CGTCTC', 'GAATTC'])"
        
        RE_sites_new = []
        for i in RE_sites:
            for k in i:
                assert k.upper() in ['A', 'T', 'C', 'G'], "RE sites must only contain A, T, C, or G (Ns and other nucleic acid codes not permitted)"
            RE_sites_new.append(i.upper())
            RE_sites_new.append(str(Bio.Seq.Seq(i).reverse_complement()).upper())

        RE_sites = RE_sites_new

        #check for u6 terminator sequence "TTTT" in pbs or protospacer
        for idx, i in pegRNA_df.iterrows():
            #check if its base or prime for filtration purposes            
            if prime==True:
                pbs_rtt = i["RTT_PBS"]
                proto = i["Protospacer"]

                if terminator_sequence in (pbs_rtt or proto):
                    contains_terminator.append(True)
                else:
                    contains_terminator.append(False)

                if sensor==True:
                    target = i['sensor_wt']
                    #finally check for restriction enzyme sites
                    res_proto = any(re_site in proto for re_site in RE_sites)
                    res_pbs_rtt = any(re_site in pbs_rtt for re_site in RE_sites)
                    res_target = any(re_site in target for re_site in RE_sites)

                    if (res_proto or res_pbs_rtt or res_target)==True:
                        RE_site.append(True)
                    else:
                        RE_site.append(False) 
                elif sensor==False:
                    #finally check for restriction enzyme sites
                    res_proto = any(re_site in proto for re_site in RE_sites)
                    res_pbs_rtt = any(re_site in pbs_rtt for re_site in RE_sites)

                    if (res_proto or res_pbs_rtt)==True:
                        RE_site.append(True)
                    else:
                        RE_site.append(False) 
                    

            elif base==True:
                proto = i["Protospacer"]
                if terminator_sequence in (proto):
                    contains_terminator.append(True)
                else:
                    contains_terminator.append(False)

                if sensor==True:
                    target = i['sensor_wt']
                    #finally check for restriction enzyme sites
                    res_proto = any(re_site in proto for re_site in RE_sites)
                    res_target = any(re_site in target for re_site in RE_sites)

                    if (res_proto or res_target)==True:
                        RE_site.append(True)
                    else:
                        RE_site.append(False) 
                elif sensor==False:
                    #finally check for restriction enzyme sites
                    res_proto = any(re_site in proto for re_site in RE_sites)
                    if (res_proto)==True:
                        RE_site.append(True)
                    else:
                        RE_site.append(False) 

        pegRNA_df['contains_polyT_terminator']=contains_terminator
        pegRNA_df['contains_RE_site']=RE_site
    
    if RE_sites == None:
        #check for u6 terminator sequence "TTTT" in pbs or protospacer
        for idx, i in pegRNA_df.iterrows():
            if prime==True:
                pbs_rtt = i["RTT_PBS"]
                proto = i["Protospacer"]

                if terminator_sequence in (pbs_rtt or proto):
                    contains_terminator.append(True)
                else:
                    contains_terminator.append(False)
            elif base==True:
                proto = i["Protospacer"]
                if terminator_sequence in (proto):
                    contains_terminator.append(True)
                else:
                    contains_terminator.append(False)
            
        pegRNA_df['contains_polyT_terminator']=contains_terminator

    return pegRNA_df

def run(input_df, input_format, chrom_dict=None, PAM = "NGG", rankby = 'PEGG2_Score', pegRNAs_per_mut = 'All',
        RTT_lengths = [5,10,15,25,30], PBS_lengths = [8,10,13,15], min_RHA_size = 1,
        RE_sites=None, polyT_threshold=4, 
        proto_size=19, context_size = 120, 
        before_proto_context=5, sensor_length=60, sensor_orientation = 'reverse-complement', sensor=True):

    """ 
    Master function for generating pegRNAs. Takes as input a dataframe containing mutations in one of the acceptable formats.
    Returns a dataframe with pegRNAs with desired design parameters.

    Parameters
    ------------
    input_df
        *type = pd.DataFrame*

        Pandas dataframe that contains the input mutations.

    input_format
        *type = str*

        Options = 'cBioPortal', 'WT_ALT', 'PrimeDesign'. For 'WT_ALT', make sure to put headers as "WT" and "ALT". 
        For "PrimeDesign", put the header as "SEQ".

    chrom_dict
        *type = dict or None*

        Dictionary generated by genome_loader() that holds the chromosome sequences.

    PAM
        *type = str*

        PAM sequence for searching. Default = "NGG". Can include any nucleic acid code (e.g. PAM = "NRCH").

    rank_by
        *type = str*

        What pegRNA parameter to rank pegRNAs by. Options = "PEGG2_Score" (weighted linear regression of different pegRNA parameters)
        or "RF_Score" (random forest predictor of pegRNA efficiency).
    
    pegRNAs_per_mut
        *type = 'All' or int*

        How many pegRNAs to produce per mutation. Default = 'All' (all possible pegRNAs with parameters). Otherwise, choose an integer value (e.g. 5).
    
    RTT_lengths
        *type = list*

        List containing RTT lengths to design pegRNAs for.

    PBS_length
        *type = list*

        List containing PBS lengths to desing pegRNAs for.

    min_RHA_size
        *type = int*

        Minimum size of the RHA (Right homology arm). Default = 1. Generally pegRNAs with smaller RHA perform poorly.

    RE_sites
        *type = list or None*

        A list containing the RE recognition sites to filter (e.g. ['CGTCTC', 'GAATTC'] for Esp3I and EcoRI). Default = None (no filtration).

    polyT_threshold
        *type = int*

        The length of the polyT sequence to classify as a terminator. Default = 4.

    proto_size
        *type = int*

        The length of the protospacer (excluding the appended G at the 5' end). Default = 19 (G+19).

    context_size
        *type = int*

        The amount of context/flanking sequence on either side of the mutation to generate. For larger variants, set this larger.
        Default = 120. e.g. AAA(A/G)AAA = context_size of 3

    before_proto_context
        *type = int*

        Default = 5. Amount of nucleotide context to put before the protospacer in the sensor

    sensor_length
        *type = int*

        Total length of the sensor in nt. Default = 60. 

    sensor_orientation
        *type = str*

        Options for sensor_orientation = 'reverse-complement' or'forward'.

    sensor
        *type = bool*

        True/False whether to include a sensor in the pegRNA design or not. 
    """


    if PAM != 'NGG':
        print('Warning: Efficiency predictions on non-NGG PAMs are un-tested')

    for i in PBS_lengths:
        assert i<=17, "Max PBS Length = 17; rerun with smaller PBS length (or use preset parameters)"

    for i in RTT_lengths:
        if i>30:
            print("Warning: RTT lengths larger than 30 nt are not reccomended due to potential low PE efficiency")
            
    ##enforce the proper formatting for INS, DEL, and INDEL
    #Note: for some reason, Clinvar encodes INS/DEL differently than cBioPortal so most INS/DEL need to be labelled as INDELs
    if input_format == 'cBioPortal':
        input_df.loc[((input_df['Variant_Type']=='INS') & (~input_df['Reference_Allele'].isin(['-', '', None]))), 'Variant_Type'] = 'INDEL'
        input_df.loc[((input_df['Variant_Type']=='DEL') & (~input_df['Tumor_Seq_Allele2'].isin(['-', '', None]))), 'Variant_Type'] = 'INDEL'

        assert chrom_dict != None, "Genome required to use cBioPortal input format. Load genome with genome_loader() and re-run."

    #and format the input df
    input_df = input_formatter(input_df, input_format, chrom_dict, context_size)

    combined_peg_dfs = []

    for i, val in input_df.iterrows():


        mut = mutation(val['wt_w_context'], val['alt_w_context'], val['left_context'], val['right_context'], val['Variant_Type'], val['REF'], val['ALT'], chrom=None, genome=None)
        mut.PAM_idx_forward, mut.PAM_idx_rc = eligible_PAM_finder(mut, PAM, max(RTT_lengths), proto_size)

        peg_df_plus = pegRNA_generator(mut, PAM, '+', proto_size, RTT_lengths, PBS_lengths)
        peg_df_minus = pegRNA_generator(mut, PAM, '-', proto_size, RTT_lengths, PBS_lengths)

        #need to translate PAM indexing for the - strand!!!!
        #convert to genome coordinates as well...
        peg_df = pd.concat((peg_df_plus, peg_df_minus))
        peg_df['mutation_idx'] = i

        combined_peg_dfs.append(peg_df)

    peg_df = pd.concat(combined_peg_dfs).set_index('mutation_idx').reset_index()

    if len(peg_df)==0:
        print('No pegRNAs generated; try increasing the size of the RTT in parameters')
        return peg_df
    
    #COMBINE WITH Input dataframe information...
    peg_df = pd.merge(input_df, peg_df, on='mutation_idx')

    #filter to exclude pegRNAs with no RHA (or whatever the input parameters are)
    peg_df = peg_df[peg_df['RHA_size']>=min_RHA_size]

    #do the pegRNA scoring
    peg_df = ontarget_score(peg_df)

    #calculate the peg_score
    peg_df = peggscore2(peg_df)

    #and calculate the RF score
    peg_df = RF_score(peg_df)

    #and then do the ranking
    uniq_muts = np.unique(peg_df['mutation_idx'])
    for mm in uniq_muts:
        subset = peg_df[peg_df['mutation_idx']==mm]
        subset = subset.sort_values(by=rankby, ascending=False)
        idxs = subset.index
        peg_df.loc[idxs, 'pegRNA_rank'] = list(range(1,len(subset)+1))

    if pegRNAs_per_mut not in ['All', 'all', ' All', ' all', 'all ', 'All ']:
        #filter
        peg_df = peg_df[peg_df['pegRNA_rank']<=pegRNAs_per_mut]

    #create sensor if desired
    if sensor == True:
        peg_df = sensor_generator(peg_df, proto_size, before_proto_context, sensor_length, sensor_orientation)
    
    #add in information about polyT sequences and RE sites
    peg_df = other_filtration(peg_df, RE_sites, polyT_threshold)

    #resetting the index and sorting by mutation
    peg_df = peg_df.sort_values(by=['mutation_idx', 'pegRNA_rank'], ascending=[True, True]).set_index('mutation_idx').reset_index()
    if 'index' in peg_df.keys():
        peg_df = peg_df.drop(columns='index')



    #restore columns to boolean
    peg_df['PAM_disrupted'] = peg_df['PAM_disrupted'].astype(bool)
    peg_df['Proto_disrupted'] = peg_df['Proto_disrupted'].astype(bool)

    #and add names...
    peg_df['pegRNA_id'] = ['pegRNA_' + str(i) for i in range(len(peg_df))]

    return peg_df


def peggscore2(df):
    """ 
    Function for calculating the PEGG2 Score. A multiple linear regression score based on pegRNA parameters.

    Parameters
    -----------
    df
        *type = pd.DataFrame*

        Dataframe containing the pegRNAs generated by run().

    """
    
    df['Proto_disrupted'] = np.array(df['Proto_disrupted'])+0
    df['PAM_disrupted'] = np.array(df['PAM_disrupted'])+0

    #also add in the PBS and RTT GC content
    #and the size of the variant
    pbs_gc = []
    rtt_gc = []
    var_size = []

    for i, val in df.iterrows():
        RTT = val['RTT']
        PBS = val['PBS']
        gc_count = RTT.count('G') + RTT.count('C')
        rtt_gc.append(gc_count/len(RTT))

        gc_count2 = PBS.count('G') + PBS.count('C')
        pbs_gc.append(gc_count2/len(PBS))

        ref = val['REF']
        alt = val['ALT']
        indelsize = abs(len(ref)-len(alt))
        var_size.append(indelsize)

    df['PBS_GC_content'] = pbs_gc
    df['RTT_GC_content'] = rtt_gc
    df['indel_size'] = var_size

    #check if there are any none value types in Ontarget score; if so set to 0
    df2 = df.copy()
    df2.loc[df2['OnTarget_Azimuth_Score'].isna(), 'OnTarget_Azimuth_Score'] = 0

    #dummy weights for PEGG Score
    factor_weights = {'Proto_disrupted':5, 
               'PAM_disrupted':3, 
               'indel_size':0, 
               'PBS_GC_content':0, 
               'RTT_GC_content':0, 
               'RTT_length':0, 
               'PBS_length':0, 
               'Distance_to_nick':-1, 
               'RHA_size':2, 
               'OnTarget_Azimuth_Score':.2}


    fw = [ 2.89996998,  1.58522483, -0.51088186,  8.78468344,  6.30826281,
       -0.33821302,  0.1058144 ,  0.30468214,  0.36783961,  0.14474713]

    factor_weights = dict(zip(factor_weights.keys(), fw))

    factor_lists = []
    for i in factor_weights.keys():
        weight = factor_weights[i]
        factor_lists.append(weight*np.asarray(df2[i]))

    pegg2_score = np.zeros(len(factor_lists[0]))
    for i in factor_lists:
        pegg2_score+=i
    
    df['PEGG2_Score']=pegg2_score

    return df

def RF_score(df):
    """ 
    Function for calculating the Random Forest score from pegRNA parameters.

    Parameters
    -----------
    df
        *type = pd.DataFrame*

        Dataframe containing the pegRNAs generated by run().
    """
    m = pickle.load(open('RF_PE.pkl', 'rb'))

    factors = ['Proto_disrupted', 'PAM_disrupted', 'indel_size', 'PBS_GC_content', 'RTT_GC_content', 'RTT_length', 'PBS_length', 'Distance_to_nick', 'RHA_size', 'OnTarget_Azimuth_Score']

    #check if there are any none value types in Ontarget score; if so set to 0
    df2 = df.copy()
    df2.loc[df2['OnTarget_Azimuth_Score'].isna(), 'OnTarget_Azimuth_Score'] = 0

    X = df2[factors]

    pred = m.predict(X)

    df['RF_Score'] = pred
    
    return df


#--------visualization tools-----------
def split_word(word):
    """
    Simple function for splitting string into component characters
    """
    return [char for char in word]

def sensor_viz(df_w_sensor, i):    
    """ 
    Function for visualizing the pegRNA aligned to the sensor. Input = dataframe with pegRNA-sensor pairs
    and the row index (i) to visualize.

    Parameters
    -----------
    df_w_sensor
        *type = pd.DataFrame*

        Dataframe containing the pegRNA-sensor pairs generated by run().

    i
        *type = int*

        Row index from the dataframe that you want to visualize.
    """
    sensor_orientation = df_w_sensor.iloc[i]['sensor_orientation']

    PAM = df_w_sensor.iloc[i]['PAM']
    PAM_len = len(PAM)

    sensor = df_w_sensor.iloc[i]['sensor_wt']
    sensor_comp = str(Bio.Seq.Seq(sensor).complement())
    split_test = split_word(sensor)
    split_test_comp = split_word(str(sensor_comp))

    if sensor_orientation == 'reverse-complement':
        #protospacer
        proto = df_w_sensor.iloc[i]['Protospacer']
        proto_start = sensor_comp.find(proto[::-1][:19])
        proto_left = ['-']*(proto_start)
        proto_right = ['-']*(len(sensor)-proto_start-len(proto))
        split_p = split_word(proto[::-1])
        split_proto = proto_left + split_p + proto_right

        #RTT_PBS
        RTT_len = df_w_sensor.iloc[i]['RTT_length']
        PBS_len =  df_w_sensor.iloc[i]['PBS_length']
        RTT_PBS = df_w_sensor.iloc[i]['RTT_PBS']
        RTT_PBS_end = proto_start+3+PBS_len
        RTT_PBS_start = RTT_PBS_end - len(RTT_PBS)
        left = ['-']*RTT_PBS_start
        right = ['-']*(len(sensor) - RTT_PBS_end)
        split_RTTPBS = split_word(RTT_PBS)
        split_RTT_PBS = left + split_RTTPBS + right

        #translation for coloring
        dict_bases = {'T':0, 'A':1, 'C':2, 'G':3, '-':4, 'X':5}
        num_translation = [dict_bases[i] for i in split_test]
        num_translation_c = [dict_bases[i] for i in split_test_comp]
        num_translation_proto = [dict_bases[i] for i in split_proto]
        num_translation_RTT = [dict_bases[i] for i in split_RTT_PBS]

        text_df = [split_proto, split_test, split_test_comp, split_RTT_PBS]
        dataFrame = [num_translation_proto, num_translation, num_translation_c,num_translation_RTT]

        fig_height = 4.5
        fig = plt.figure(figsize=(20,fig_height))
        myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white')
        #if var_type=='DEL':
        #    myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white', 'black')

        cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

        ax = sns.heatmap(dataFrame, annot=text_df, fmt="", linewidth=0, cmap=cmap, cbar=False,linewidths=0, linecolor='lightgray', xticklabels=False)

        ax.add_patch(patches.Rectangle((proto_start-PAM_len, 2), PAM_len, 1, fill=False, edgecolor='tab:green', lw=3, label='PAM')) #protospacer loca)) #protospacer location

        #protospacer
        ax.add_patch(patches.Rectangle((proto_start, 0), len(proto), 1, fill=False, edgecolor='tab:blue', lw=3, label='Protospacer')) #protospacer location

        #PBS
        ax.add_patch(patches.Rectangle((proto_start+3, 3), PBS_len, 1, fill=False, edgecolor='tab:purple', lw=3, label='PBS')) #protospacer location

        #RTT
        ax.add_patch(patches.Rectangle((proto_start+3-RTT_len, 3), RTT_len, 1, fill=False, edgecolor='yellow', lw=3, label='RTT')) #protospacer location
        ax.add_patch(Polygon([(proto_start+3-RTT_len, 3), (proto_start+3-RTT_len-2, 3.5), (proto_start+3-RTT_len, 4)],facecolor='black',edgecolor='yellow', lw=3))

        #finally, mark the location of the sequence to be mutated
        RHA = df_w_sensor.iloc[i]['RHA_size']
        vt = df_w_sensor.iloc[i]['Variant_Type']
        mut_start = proto_start+3-RTT_len + RHA

        if vt == 'INS':
            #mut_size = 2
            mut_size = len(df_w_sensor.iloc[i]['ALT'])

            ax.add_patch(patches.Rectangle((mut_start+mut_size, 1), 0, 2, fill=False, edgecolor='tab:red', lw=3, label='Mutated Bases')) #protospacer location

        elif vt=='DEL':
            mut_size = len(df_w_sensor.iloc[i]['REF'])

            ax.add_patch(patches.Rectangle((mut_start-mut_size, 1), mut_size, 2, fill=False, edgecolor='tab:red', lw=3, label='Mutated Bases')) #protospacer location


        else:
            mut_size = len(df_w_sensor.iloc[i]['REF'])

            ax.add_patch(patches.Rectangle((mut_start, 1), mut_size, 2, fill=False, edgecolor='tab:red', lw=3, label='Mutated Bases')) #protospacer location

        #also label the RHA

        ax.add_patch(patches.Rectangle((proto_start+3-RTT_len, 3), RHA, 1, fill=False, edgecolor='tab:orange', lw=3, label='RHA')) #protospacer location

    elif sensor_orientation == 'forward':

        #protospacer
        proto = df_w_sensor.iloc[i]['Protospacer']
        proto_start = sensor.find(proto[1:])-1 #account for G+19 or G+20
        proto_left = ['-']*(proto_start)
        proto_right = ['-']*(len(sensor)-proto_start-len(proto))
        split_p = split_word(proto)
        split_proto = proto_left + split_p + proto_right

        #RTT_PBS
        RTT_len = df_w_sensor.iloc[i]['RTT_length']
        PBS_len =  df_w_sensor.iloc[i]['PBS_length']
        RTT_PBS = str(df_w_sensor.iloc[i]['RTT_PBS'])[::-1] #[::-1] #reverse it for 3' to

        RTT_PBS_end = proto_start+3+PBS_len
        #RTT_PBS_start = RTT_PBS_end - len(RTT_PBS)

        RTT_PBS_start = proto_start + len(proto) -3 - PBS_len
        RTT_PBS_end = proto_start + len(proto) -3 + RTT_len

        left = ['-']*RTT_PBS_start
        right = ['-']*(len(sensor) - RTT_PBS_end)
        split_RTTPBS = split_word(RTT_PBS)
        split_RTT_PBS = left + split_RTTPBS + right

        #translation for coloring
        dict_bases = {'T':0, 'A':1, 'C':2, 'G':3, '-':4, 'X':5}
        num_translation = [dict_bases[i] for i in split_test]
        num_translation_c = [dict_bases[i] for i in split_test_comp]
        num_translation_proto = [dict_bases[i] for i in split_proto]
        num_translation_RTT = [dict_bases[i] for i in split_RTT_PBS]

        text_df = [split_RTT_PBS, split_test, split_test_comp, split_proto]
        dataFrame = [num_translation_RTT, num_translation, num_translation_c, num_translation_proto,]

        fig_height = 4.5
        fig = plt.figure(figsize=(20,fig_height))
        myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white')
        #if var_type=='DEL':
        #    myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white', 'black')

        cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

        ax = sns.heatmap(dataFrame, annot=text_df, fmt="", linewidth=0, cmap=cmap, cbar=False,linewidths=0, linecolor='lightgray', xticklabels=False)

        #PAM
        ax.add_patch(patches.Rectangle((proto_start+len(proto), 1), PAM_len, 1, fill=False, edgecolor='tab:green', lw=3, label='PAM')) #protospacer loca)) #protospacer location

        #protospacer
        ax.add_patch(patches.Rectangle((proto_start, 3), len(proto), 1, fill=False, edgecolor='tab:blue', lw=3, label='Protospacer')) #protospacer location

        #PBS
        ax.add_patch(patches.Rectangle((RTT_PBS_start, 0), PBS_len, 1, fill=False, edgecolor='tab:purple', lw=3, label='PBS')) #protospacer location

        #RTT
        ax.add_patch(patches.Rectangle((RTT_PBS_start + PBS_len, 0), RTT_len, 1, fill=False, edgecolor='yellow', lw=3, label='RTT')) #protospacer location
        ax.add_patch(Polygon([(RTT_PBS_start + PBS_len+RTT_len, 0), (RTT_PBS_start + PBS_len+RTT_len+2, 0.5), (RTT_PBS_start + PBS_len+RTT_len, 1)],facecolor='black',edgecolor='yellow', lw=3))

        #finally, mark the location of the sequence to be mutated
        RHA = df_w_sensor.iloc[i]['RHA_size']
        vt = df_w_sensor.iloc[i]['Variant_Type']
        mut_start = RTT_PBS_start + RTT_len + PBS_len - RHA -1

        if vt == 'INS':
            #mut_size = 2
            mut_size = len(df_w_sensor.iloc[i]['ALT'])

            ax.add_patch(patches.Rectangle((mut_start-mut_size+1, 1), 0, 2, fill=False, edgecolor='tab:red', lw=3, label='Mutated Bases')) #protospacer location

        elif vt=='DEL':
            mut_size = len(df_w_sensor.iloc[i]['REF'])

            ax.add_patch(patches.Rectangle((mut_start-mut_size, 1), mut_size, 2, fill=False, edgecolor='tab:red', lw=3, label='Mutated Bases')) #protospacer location


        else:
            mut_size = len(df_w_sensor.iloc[i]['REF'])

            ax.add_patch(patches.Rectangle((mut_start-mut_size+1, 1), mut_size, 2, fill=False, edgecolor='tab:red', lw=3, label='Mutated Bases')) #protospacer location

        #also label the RHA

        ax.add_patch(patches.Rectangle((mut_start+1, 0), RHA, 1, fill=False, edgecolor='tab:orange', lw=3, label='RHA')) #protospacer location

    ax.set_yticklabels(["3'", "5'", "3'", "5'"], rotation=0, fontsize=14)

    ax.legend(bbox_to_anchor=(1.16,.78), fontsize=14)

    ref_allele = df_w_sensor.iloc[i]['REF']
    mut_allele = df_w_sensor.iloc[i]['ALT']

    ax.set_title('RTT length: '+str(RTT_len)+ ' nt, PBS length: ' + str(PBS_len)+' nt | ' + vt + ': '+ref_allele + '>' + mut_allele, fontsize=16)


#----Cloning---------
def goldengate_oligos(peg_df):
    """ 
    Currently a place-holder. Going to implement automated Golden Gate oligo generator (so scaffold doesn't need to be synthesized).
    """
    return None

def prime_oligo_generator(peg_df, epeg=True, epeg_motif = 'tevopreQ1',
                    five_prime_adapter = 'AGCGTACACGTCTCACACC',three_prime_adapter = 'GAATTCTAGATCCGGTCGTCAAC',
                gRNA_scaff = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'):

    """
    A tool for automatically generating oligos from the output of run().
    Returns input dataframe with new columns containing the pegRNA oligo.

    Parameters
    ----------
    peg_df
        *type = pd.DataFrame*
        
        A dataframe containing the pegRNAs for the selected input mutations. Generated by run() or pegRNA_generator()

    epeg
        *type = bool*

        True/False whether to include an epegRNA motif at the end of the 3' extension.

    epeg_motif
        *type = str*

        Which epeg motif to include at end of 3' extension. Default = "tevopreQ1"; Other option = "mpknot" (both from Chen et al.).
        Can also input a custom motif sequence (just put in the sequence in 5' to 3' orientation)!
    
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
    
    """
        
        
    u6_term = 'TTTTTTT'
    tevopreQ1 = 'CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA'
    mpknot = "GGGTCAGGAGCCCCCCCCCTGAACCCAGGATAACCCTCAAAGTCGGGGGGCAACCC"

    if epeg==True:
        if epeg_motif == 'tevopreQ1':
            epeg_motif = tevopreQ1

        elif epeg_motif == 'mpknot':
            epeg_motif = mpknot

        else:
            epeg_motif=epeg_motif
            for i in epeg_motif:
                assert i in ['A','T','C','G','a','t','c','g'], "Choose a valid epeg_motif. Your custom epeg_motif contains non-DNA bases."

    elif epeg==False:
        epeg_motif = ''

    peg_oligos = []
    #check if sensor is in the peg_df
    if "sensor_wt" in peg_df.keys():
        for i, val in peg_df.iterrows():
            proto = val['Protospacer']
            extension = val['RTT_PBS']
            sensor = val["sensor_wt"]
            if sensor==None: #there are some sensor errors where it is =None
                sensor = ''

            pegRNA_full = five_prime_adapter + proto + gRNA_scaff + extension + epeg_motif + u6_term + sensor + three_prime_adapter
            peg_oligos.append(pegRNA_full)

    #if it isn't:
    else:
        for i, val in peg_df.iterrows():
            proto = val['Protospacer']
            extension = val['RTT_PBS']

            pegRNA_full = five_prime_adapter + proto + gRNA_scaff + extension + epeg_motif + u6_term + three_prime_adapter
            peg_oligos.append(pegRNA_full)


    peg_df['pegRNA_oligo'] = peg_oligos
    
    return peg_df
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import Bio.Seq
import gzip
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.pairwise2 import format_alignment
import warnings
import regex as re
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
    Takes in filepath of human genome (GrCH37 or GrCh38) and returns records and index_list for PEGG parsing.
    
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

def df_formatter(df, chrom_dict, context_size = 120, 
                 cols_to_save = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Consequence', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2',  'HGVSc', 'HGVSp_Short']):

    """ 
    Takes in variants (in cBioPortal format!)
    and outputs dataframe with WT and ALT oligos with designated context_size
    context_size = the amount of nt on either side of the variant e.g. AAA(A/G)AAA = context_size of 3
    cols_to_save = information to save about the mutations...(can modify beyond defaults...)
    """

    wt_w_context = []
    alt_w_context = []
    left_context_list = []
    right_context_list = []
    ref_allele = []
    alt_allele = []

    seq_start = []
    seq_end = []

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

        if vt in ['SNP', 'ONP', 'DNP']:
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

        wt_w_context.append(str(wt_seq))
        alt_w_context.append(str(alt_seq))
        left_context_list.append(str(left_context))
        right_context_list.append(str(right_context))
        ref_allele.append(str(ref))
        alt_allele.append(str(alt))

        start = s-context_size
        end = e+context_size

        seq_start.append(start)
        seq_end.append(end)

        assert str(chr_seq[start-1:end])==str(wt_seq), print(chr_seq[start-1:end] + '\n' + str(wt_seq))
                                                            

    df_new = df[cols_to_save]

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
    This program will break for more complex mutations (e.g. INDELs)
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
    takes as input sequence in prime design format 
    e.g. AATTCCG(G/C)AATTCGCT
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

def PAM_finder(seq, PAM):
    """ 
    Finds indeces of PAM sequences
    Formatting: e.g. "NGG"
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
    Determines eligible PAM sequences for creating pegRNAs...

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

def sensor_generator(df, before_proto_context=5, sensor_length=60, sensor_orientation = 'reverse-complement'):
    """ 
    Generates sensor sequence for quantification of pegRNA editing outcomes
    Automatically puts sensor in reverse complement orientation with respect to protospacer
    This is highly reccomended to reduce recombination during cloning/library preparation

    options for sensor_orienation = 'reverse_comp'; 'forward'

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

        s1 = pam_start-20-before_proto_context
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

    return df


def input_formatter(input_df, input_format, chrom_dict, context_size):

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

    elif input_format == 'ClinVar':
        return "ClinVar format input in development"
        #ADD FUNCTIONS...
    
    return input_df


def run(input_df, input_format, chrom_dict, PAM = "NGG", RTT_lengths = [5,10,15,25,30], PBS_lengths = [8,10,13,15], 
        proto_size=19, context_size = 120, 
        before_proto_context=5, sensor_length=60, sensor_orientation = 'reverse-complement', sensor=True):

    for i in PBS_lengths:
        assert i<=17, "Max PBS Length = 17; rerun with smaller PBS length (or use preset parameters)"

    for i in RTT_lengths:
        if i>30:
            print("Warning: RTT lengths larger than 30 nt are not reccomended due to potential low PE efficiency")
            
    #format the input df
    #NEED TO UPDATE TO INCLUDE CLINVAR!!!
    #also alter the "cols_to_save" parameter to avoid errors here...
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

    #COMBINE WITH Input dataframe information...
    peg_df = pd.merge(input_df, peg_df, on='mutation_idx')

    #do the pegRNA scoring
    #first add in the on-target score...
    #MAYBE LIMIT THIS TO ONLY NGG PAMS???
    peg_df = ontarget_score(peg_df)

    #then calculate the PEGG Score and filter out the pegRNAs...
    #have an option for "All" possible pegRNAs with given parameters

    #create sensor if desired
    if sensor == True:
        peg_df = sensor_generator(peg_df, before_proto_context, sensor_length, sensor_orientation)
    
    return peg_df


def peggscore2(df):
    
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
        factor_lists.append(weight*np.asarray(df[i]))

    pegg2_score = np.zeros(len(factor_lists[0]))
    for i in factor_lists:
        pegg2_score+=i
    
    df['PEGG2_Score']=pegg2_score

    return df
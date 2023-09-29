from pegg.prime import *
import numpy as np
import pandas as pd
import Bio.Seq
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import matplotlib.patches as patches
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

def eligible_PAM_finder_base(mut, PAM, proto_size=19):
    """ 
    Determines eligible PAM sequences for creating gRNAs.
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

    proto_size
        *type = int*

        Size of protospacer being used. Default = 19 (G+19).
    """

    forward = pegg2.PAM_finder(mut.wt_forward, PAM)
    reverse = pegg2.PAM_finder(mut.wt_rc, PAM)

    var_type = mut.variant_type

    #first do the forward sequence
    left_len_F = len(mut.left_seq)
    left_len_R = len(mut.left_seq_rc)


    eligible_PAM_start_F = left_len_F + 1  #
    eligible_PAM_end_F = left_len_F + 1 + 19 
    
    eligible_PAM_start_R = left_len_R + 1
    eligible_PAM_end_R = left_len_R + 1 + 19 

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


def gRNA_generator(mut, PAM, orientation, proto_size, ideal_edit_window = [4,8]):
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

    ideal_edit_window
        *type = list*

        Ideal editing window for the editor being used. Default = [4,8]. Labels mutations that fall in this window for 
        future filtration if desired.
    """

    #------function----
    protospacer_seqs = []
    protospacer_wide_seqs = []
    PAM_start_list = []
    PAM_end_list = []
    PAM_list = []
    proto_location = []
    ideal_edit = []

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

        proto_loc = PAM_start - left_len
        #and then convert to 1 to 20 indexing...
        proto_loc = 21-proto_loc
        proto_location.append(proto_loc)
        if (proto_loc>=ideal_edit_window[0]) and (proto_loc<=ideal_edit_window[0]):
            ideal_edit.append(True)
        else:
            ideal_edit.append(False)

    #pull out the PAM sequence and protospacer as well
        PAM_sequence = seq_F[PAM_start:PAM_end]
        protospacer =  'G' + seq_F[PAM_start - proto_size:PAM_start]

        #extract wide protospacer for Rule set 2 prediction...
        protospacer_wide = seq_F[PAM_start - 20 - 4 : PAM_start + 3 + 3]

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

    cols = [PAM_start_list, PAM_end_list, PAM_list, orientation, protospacer_wide_seqs, protospacer_seqs, proto_location, ideal_edit]
    col_labels = ["PAM_start", "PAM_end", "PAM", "PAM_strand","Protospacer_30", "Protospacer", "Protospacer_Location", "Ideal_Edit_Window"]
    dtypes = ['int', 'int', 'str', 'str', 'str', 'str', 'int', 'bool']

    dtype_dict = dict(zip(col_labels, dtypes))

    peg_df = pd.DataFrame(dict(zip(col_labels, cols))).dropna().reset_index().drop(columns='index').astype(dtype_dict)

    return peg_df



def run_base(input_df, input_format, chrom_dict=None, PAM = "NGG", filtration = "ABE+CBE", ideal_edit_window = [4,8], auto_SNP_filter = True, 
        proto_size=19, context_size = 120, 
        RE_sites=None, polyT_threshold=4, 
        before_proto_context=5, sensor_length=40, sensor_orientation = 'reverse-complement', sensor=True):
            
    """ 
    Master function for generating base editing gRNAs. Takes as input a dataframe containing mutations in one of the acceptable formats.
    Returns a dataframe with gRNAs with desired design parameters.

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

    filtration
        *type = str or list*

        Filters the mutation input list to only include the desired SNPs. Options = "No filter", "ABE", "CBE", "ABE+CBE",
        or a list containing the desired SNPs to model (e.g. ['C>A', 'T>C']).

    ideal_edit_window
        *type = list*

        Ideal editing window for the editor being used. Default = [4,8]. Labels mutations that fall in this window for 
        future filtration if desired.

    auto_SNP_filter
        *type = bool*

        True/False for whether to filter mutant input to exclude mutations that are NOT SNPs (and thus not BE amenable).
    
    proto_size
        *type = int*

        The length of the protospacer (excluding the appended G at the 5' end). Default = 19 (G+19).

    context_size
        *type = int*

        The amount of context/flanking sequence on either side of the mutation to generate. For larger variants, set this larger.
        Default = 120. e.g. AAA(A/G)AAA = context_size of 3

    RE_sites
        *type = list or None*

        A list containing the RE recognition sites to filter (e.g. ['CGTCTC', 'GAATTC'] for Esp3I and EcoRI). Default = None (no filtration).

    polyT_threshold
        *type = int*

        The length of the polyT sequence to classify as a terminator. Default = 4.

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

    input_df = pegg2.input_formatter(input_df, input_format, chrom_dict, context_size)

    vts = np.unique(input_df['Variant_Type'])

    #make sure it's a SNP
    if auto_SNP_filter==False:
        not_allowed = ['ONP', 'DEL', 'INS', 'INDEL']
        for i in not_allowed:
            assert i not in vts, "Only SNPs can be modeled with base editing; remove non-SNPs from mutation list or select auto_SNP_filter=True"

    elif auto_SNP_filter==True:
        input_df = input_df[input_df['Variant_Type']=='SNP']
        not_allowed = ['ONP', 'DEL', 'INS', 'INDEL']
        for i in not_allowed:
            if i in vts:
                print("Only SNPs can be modeled with base editing; Non-SNPs automatically removed from list")

    combined_peg_dfs = []

    for i, val in input_df.iterrows():


        mut = pegg2.mutation(val['wt_w_context'], val['alt_w_context'], val['left_context'], val['right_context'], val['Variant_Type'], val['REF'], val['ALT'], chrom=None, genome=None)
        mut.PAM_idx_forward, mut.PAM_idx_rc = eligible_PAM_finder_base(mut, PAM, proto_size)

        peg_df_plus = gRNA_generator(mut, PAM, '+', proto_size, ideal_edit_window)
        peg_df_minus = gRNA_generator(mut, PAM, '-', proto_size, ideal_edit_window)

        #need to translate PAM indexing for the - strand!!!!
        #convert to genome coordinates as well...
        peg_df = pd.concat((peg_df_plus, peg_df_minus))
        peg_df['mutation_idx'] = i

        combined_peg_dfs.append(peg_df)

    peg_df = pd.concat(combined_peg_dfs).set_index('mutation_idx').reset_index()

    #COMBINE WITH Input dataframe information...
    peg_df = pd.merge(input_df, peg_df, on='mutation_idx')

    #and then do filtration...
    if filtration != 'No filter':
        if filtration == 'CBE':
            p1 = peg_df[(peg_df['REF']=='C') & (peg_df['ALT']=='T') & (peg_df['PAM_strand']=='+')]
            p2 = peg_df[(peg_df['REF']=='G') & (peg_df['ALT']=='A') & (peg_df['PAM_strand']=='-')]
            peg_df = pd.concat((p1, p2))
            peg_df['Editor'] = 'CBE'

        elif filtration == 'ABE':
            p1 = peg_df[(peg_df['REF']=='A') & (peg_df['ALT']=='G') & (peg_df['PAM_strand']=='+')]
            p2 = peg_df[(peg_df['REF']=='T') & (peg_df['ALT']=='C') & (peg_df['PAM_strand']=='-')]
            peg_df = pd.concat((p1, p2))
            peg_df['Editor'] = 'ABE'

        elif filtration in ['ABE+CBE', 'CBE+ABE', 'ABE + CBE', 'CBE + ABE']:
            #filter for only ABE/CBE eligible variants and PAM strand...
            p1 = peg_df[(peg_df['REF']=='C') & (peg_df['ALT']=='T') & (peg_df['PAM_strand']=='+')]
            p2 = peg_df[(peg_df['REF']=='G') & (peg_df['ALT']=='A') & (peg_df['PAM_strand']=='-')]
            peg_df1 = pd.concat((p1, p2))
            peg_df1['Editor'] = 'CBE'

            p1 = peg_df[(peg_df['REF']=='A') & (peg_df['ALT']=='G') & (peg_df['PAM_strand']=='+')]
            p2 = peg_df[(peg_df['REF']=='T') & (peg_df['ALT']=='C') & (peg_df['PAM_strand']=='-')]
            peg_df2 = pd.concat((p1, p2))
            peg_df2['Editor'] = 'ABE'

            peg_df = pd.concat((peg_df1, peg_df2))

        elif type(filtration)==list:

            p_holder = []
            for k in filtration:
                spl = k.split('>')
                r = spl[0]
                r_rc = str(Bio.Seq.Seq(r).complement())
                a = spl[1] 
                a_rc = str(Bio.Seq.Seq(a).complement())

                p1 = peg_df[(peg_df['REF']==r) & (peg_df['ALT']==a) & (peg_df['PAM_strand']=='+')]
                p2 = peg_df[(peg_df['REF']==r_rc) & (peg_df['ALT']==a_rc) & (peg_df['PAM_strand']=='-')]
                peg_df2 = pd.concat((p1, p2))
                peg_df2['Editor'] = k
                p_holder.append(peg_df2)

            peg_df = pd.concat(p_holder)

        else:
            assert 1==2, "Choose valid filtration format ('ABE', 'CBE', 'ABE+CBE', or list format (['C>G', ...]))"


    #resetting the index
    peg_df = peg_df.sort_values(by='mutation_idx', ascending=True).set_index('mutation_idx').reset_index()
    if 'index' in peg_df.keys():
        peg_df = peg_df.drop(columns='index')

    #and add names...
    peg_df['gRNA_id'] = ['gRNA_' + str(i) for i in range(len(peg_df))]

    #do the pegRNA scoring
    #first add in the on-target score...
    #MAYBE LIMIT THIS TO ONLY NGG PAMS???
    peg_df = pegg2.ontarget_score(peg_df)

    #then calculate the PEGG Score and filter out the pegRNAs...
    #have an option for "All" possible pegRNAs with given parameters

    #create sensor if desired
    if sensor == True:
        peg_df = sensor_generator_base(peg_df, proto_size, before_proto_context, sensor_length, sensor_orientation)

    #add in information about polyT sequences and RE sites
    peg_df = pegg2.other_filtration(peg_df, RE_sites, polyT_threshold)
    
    return peg_df


def sensor_generator_base(df, proto_size, before_proto_context=5, sensor_length=60, sensor_orientation = 'reverse-complement'):
    """ 
    Generates sensor sequence for quantification of gRNA editing outcomes.
    Automatically puts sensor in reverse complement orientation with respect to protospacer.
    This is highly reccomended to reduce recombination during cloning/library preparation.

    Parameters
    -----------
    df
        *type = pd.DataFrame*

        Dataframe containing gRNAs.

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
        #RTT_PBS = val['RTT_PBS']
        #RHA = val['RHA_size']

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

                   
                        sensor_wt_seqs.append(sensor_wt)
                        sensor_alt_seqs.append(sensor_alt)
                        sensor_error.append("No Error")
                    
                    if sensor_orientation == 'forward':
                    
                        sensor_wt_seqs.append(sensor_wt)
                        sensor_alt_seqs.append(sensor_alt)
                        sensor_error.append("No Error")


    df['sensor_wt'] = sensor_wt_seqs
    df['sensor_alt'] = sensor_alt_seqs
    df['sensor_orientation'] = sensor_orientation
    df['sensor_error'] = sensor_error

    return df

#--------Visualization tools----------
def split_word(word):
    """
    Simple function for splitting string into component characters
    """
    return [char for char in word]

def sensor_viz_base(df_w_sensor, i):    
    """ 
    Function for visualizing the gRNA aligned to the sensor. Input = dataframe with pegRNA-sensor pairs
    and the row index (i) to visualize.

    Parameters
    -----------
    df_w_sensor
        *type = pd.DataFrame*

        Dataframe containing the gRNA-sensor pairs generated by run().

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

    proto_loc = df_w_sensor.iloc[i]['Protospacer_Location']

    if sensor_orientation == 'reverse-complement':
        #protospacer
        proto = df_w_sensor.iloc[i]['Protospacer']
        proto_start = sensor_comp.find(proto[::-1][:19])
        proto_left = ['-']*(proto_start)
        proto_right = ['-']*(len(sensor)-proto_start-len(proto))
        split_p = split_word(proto[::-1])
        split_proto = proto_left + split_p + proto_right

        #translation for coloring
        dict_bases = {'T':0, 'A':1, 'C':2, 'G':3, '-':4, 'X':5}
        num_translation = [dict_bases[i] for i in split_test]
        num_translation_c = [dict_bases[i] for i in split_test_comp]
        num_translation_proto = [dict_bases[i] for i in split_proto]

        text_df = [split_proto, split_test, split_test_comp]
        dataFrame = [num_translation_proto, num_translation, num_translation_c]

        fig_height = 4.5
        fig = plt.figure(figsize=(20,fig_height))
        myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white')
        #if var_type=='DEL':
        #    myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white', 'black')

        cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

        ax = sns.heatmap(dataFrame, annot=text_df, fmt="", linewidth=0, cmap=cmap, cbar=False,linewidths=0, linecolor='lightgray', xticklabels=False, yticklabels=True)
        ax.set_yticklabels(["3'", "5'", "3'"], rotation=0, fontsize=14)


        ax.add_patch(patches.Rectangle((proto_start-PAM_len, 2), PAM_len, 1, fill=False, edgecolor='tab:green', lw=3, label='PAM')) #protospacer loca)) #protospacer location

        #protospacer
        ax.add_patch(patches.Rectangle((proto_start, 0), len(proto), 1, fill=False, edgecolor='tab:blue', lw=3, label='Protospacer')) #protospacer location

        #mutant location
        mut_start = proto_start+(20-proto_loc)
        print(mut_start)
        ax.add_patch(patches.Rectangle((mut_start, 2), 1, 1, fill=False, edgecolor='tab:red', lw=3, label='Mutated Base')) #protospacer location

    elif sensor_orientation == 'forward':

        #protospacer
        proto = df_w_sensor.iloc[i]['Protospacer']
        proto_start = sensor.find(proto[1:])-1 #account for G+19 or G+20
        proto_left = ['-']*(proto_start)
        proto_right = ['-']*(len(sensor)-proto_start-len(proto))
        split_p = split_word(proto)
        split_proto = proto_left + split_p + proto_right

        
        #translation for coloring
        dict_bases = {'T':0, 'A':1, 'C':2, 'G':3, '-':4, 'X':5}
        num_translation = [dict_bases[i] for i in split_test]
        num_translation_c = [dict_bases[i] for i in split_test_comp]
        num_translation_proto = [dict_bases[i] for i in split_proto]

        text_df = [split_test, split_test_comp, split_proto]
        dataFrame = [num_translation, num_translation_c, num_translation_proto]

        fig_height = 4.5
        fig = plt.figure(figsize=(20,fig_height))
        myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white')
        #if var_type=='DEL':
        #    myColors = ('plum', 'skyblue', 'navajowhite', 'lightgoldenrodyellow', 'white', 'black')

        cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

        ax = sns.heatmap(dataFrame, annot=text_df, fmt="", linewidth=0, cmap=cmap, cbar=False,linewidths=0, linecolor='lightgray', xticklabels=False)

        ax.set_yticklabels(["5'", "3'", "5'"], rotation=0, fontsize=14)

        #PAM
        ax.add_patch(patches.Rectangle((proto_start+len(proto), 0), PAM_len, 1, fill=False, edgecolor='tab:green', lw=3, label='PAM')) #protospacer loca)) #protospacer location

        #protospacer
        ax.add_patch(patches.Rectangle((proto_start, 2), len(proto), 1, fill=False, edgecolor='tab:blue', lw=3, label='Protospacer')) #protospacer location

        #mutation_location
        mut_start = proto_loc+proto_start-1
        ax.add_patch(patches.Rectangle((mut_start, 0), 1, 1, fill=False, edgecolor='tab:red', lw=3, label='Mutated Base')) #protospacer location


    ax.legend(bbox_to_anchor=(1.16,.78), fontsize=14)

    ref_allele = df_w_sensor.iloc[i]['REF']
    mut_allele = df_w_sensor.iloc[i]['ALT']
    strand = df_w_sensor.iloc[i]['PAM_strand']

    ax.set_title(f'PAM strand = {strand} | ' + ref_allele + '>' + mut_allele + f" | Protospacer Location = +{proto_loc}", fontsize=16)


#------oligo generation tool
def base_oligo_generator(peg_df, five_prime_adapter = 'AGCGTACACGTCTCACACC',three_prime_adapter = 'GAATTCTAGATCCGGTCGTCAAC',
                gRNA_scaff = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'):

    """
    A tool for automatically generating oligos from the output of run().
    Returns input dataframe with new columns containing gRNA oligo with or without sensor.

    Parameters
    ----------
    peg_df
        *type = pd.DataFrame*
        
        A dataframe containing the gRNAs for the selected input mutations. Generated by run_base() or gRNA_generator()

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

    base_oligos = []
    #check if sensor is in the peg_df
    if "sensor_wt" in peg_df.keys():
        for i, val in peg_df.iterrows():
            proto = val['Protospacer']
            extension = val['RTT_PBS']
            sensor = val["sensor_wt"]
            if sensor==None: #there are some sensor errors where it is =None
                sensor = ''

            gRNA_full = five_prime_adapter + proto + gRNA_scaff + u6_term + sensor + three_prime_adapter
            base_oligos.append(pegRNA_full)

    #if it isn't:
    else:
        for i, val in peg_df.iterrows():
            proto = val['Protospacer']
            extension = val['RTT_PBS']

            gRNA_full = five_prime_adapter + proto + gRNA_scaff + u6_term + three_prime_adapter
            base_oligos.append(pegRNA_full)


    peg_df['gRNA_oligo'] = base_oligos
    
    return peg_df
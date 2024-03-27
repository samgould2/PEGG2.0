
import numpy as np
import os
import pandas
import sklearn
import pickle
import time
import Bio.SeqUtils as SeqUtil
import Bio.Seq as Seq
import sys
import Bio.SeqUtils.MeltingTemp as Tm
import itertools
from math import exp       
from re import findall 
#import importlib
from importlib.resources import files


##code taken from: https://github.com/maximilianh/crisporWebsite
##Modified to run within this package...

def calcAziScore(seqs):
    " the official implementation of the Doench2016 (aka Fusi) score from Microsoft "
    res = []
    for seq in seqs:
        if "N" in seq:
            res.append(-1) # can't do Ns
            continue

        pam = seq[25:27]
        # pam_audit = do not check for NGG PAM
        seq = seq.upper()
        score = predict(np.array([seq]), None, None, pam_audit=False)[0]
        res.append(int(round(100*score)))
    return res

#------function from model_comparison.py------------

def predict(seq, aa_cut=-1, percent_peptide=-1, model=None, model_file=None, pam_audit=True, length_audit=False, learn_options_override=None):
    """
    if pam_audit==False, then it will not check for GG in the expected position
    this is useful if predicting on PAM mismatches, such as with off-target
    """
    # assert not (model is None and model_file is None), "you have to specify either a model or a model_file"
    assert isinstance(seq, (np.ndarray)), "Please ensure seq is a numpy array"
    assert len(seq[0]) > 0, "Make sure that seq is not empty"
    assert isinstance(seq[0], str), "Please ensure input sequences are in string format, i.e. 'AGAG' rather than ['A' 'G' 'A' 'G'] or alternate representations"

    if aa_cut is not None:
        assert len(aa_cut) > 0, "Make sure that aa_cut is not empty"
        assert isinstance(aa_cut, (np.ndarray)), "Please ensure aa_cut is a numpy array"
        assert np.all(np.isreal(aa_cut)), "amino-acid cut position needs to be a real number"

    if percent_peptide is not None:
        assert len(percent_peptide) > 0, "Make sure that percent_peptide is not empty"
        assert isinstance(percent_peptide, (np.ndarray)), "Please ensure percent_peptide is a numpy array"
        assert np.all(np.isreal(percent_peptide)), "percent_peptide needs to be a real number"


    if model_file is None:

        ##NEED TO CHANGE THESE WITH IMPORTLIB!!!!!!

        #azimuth_saved_model_dir = os.path.join(os.path.dirname(azimuth.__file__), 'saved_models')

        if np.any(percent_peptide == -1) or (percent_peptide is None and aa_cut is None):
            #print("No model file specified, using V3_model_nopos")
            model_name = 'V3_model_nopos.pickle'
            #model_file = importlib.resources.files(__package__).joinpath(model_name)
            model_file = files(__package__).joinpath(model_name)
            
        else:
            #print("No model file specified, using V3_model_full")
            model_name = 'V3_model_full.pickle'
            #model_file = importlib.resources.files(__package__).joinpath(model_name)
            model_file = files(__package__).joinpath(model_name)

        #model_file = os.path.join(azimuth_saved_model_dir, model_name)

    if model is None:
        with open(model_file, 'rb') as f:
            model, learn_options = pickle.load(f)
    else:
        model, learn_options = model
        
    learn_options["V"] = 2

    learn_options = override_learn_options(learn_options_override, learn_options)

    # Y, feature_sets, target_genes, learn_options, num_proc = setup(test=False, order=2, learn_options=learn_options, data_file=test_filename)
    # inputs, dim, dimsum, feature_names = pd.concatenate_feature_sets(feature_sets)

    Xdf = pandas.DataFrame(columns=['30mer', 'Strand'], data=list(zip(seq, ['NA' for x in range(len(seq))])))

    if np.all(percent_peptide != -1) and (percent_peptide is not None and aa_cut is not None):
        gene_position = pandas.DataFrame(columns=['Percent Peptide', 'Amino Acid Cut position'], data=list(zip(percent_peptide, aa_cut)))
    else:
        gene_position = pandas.DataFrame(columns=['Percent Peptide', 'Amino Acid Cut position'], data=list(zip(np.ones(seq.shape[0])*-1, np.ones(seq.shape[0])*-1)))

    feature_sets = featurize_data(Xdf, learn_options, pandas.DataFrame(), gene_position, pam_audit=pam_audit, length_audit=length_audit)
    inputs, dim, dimsum, feature_names = concatenate_feature_sets(feature_sets)
    
    #print "CRISPR"
    #pandas.DataFrame(inputs).to_csv("CRISPR.inputs.test.csv")
    #import ipdb; ipdb.set_trace()

    # call to scikit-learn, returns a vector of predicted values    
    preds = model.predict(inputs)

    # also check that predictions are not 0/1 from a classifier.predict() (instead of predict_proba() or decision_function())
    unique_preds = np.unique(preds)
    ok = False
    for pr in preds:
        if pr not in [0,1]:
            ok = True
    assert ok, "model returned only 0s and 1s"
    return preds

def override_learn_options(learn_options_override, learn_options):
    """
    override all keys seen in learn_options_override to alter learn_options
    """
    if learn_options_override is not None:
        for k in list(learn_options_override.keys()):
            learn_options[k] = learn_options_override[k]
    return learn_options

def concatenate_feature_sets(feature_sets, keys=None):
    '''
    Given a dictionary of sets of features, each in a Pandas.DataFrame,
    concatenate them together to form one big np.array, and get the dimension
    of each set
    Returns: inputs, dim
    '''
    assert feature_sets != {}, "no feature sets present"
    if keys is None:
        keys = list(feature_sets.keys())

    F = feature_sets[keys[0]].shape[0]
    for set in list(feature_sets.keys()):
        F2 = feature_sets[set].shape[0]
        assert F == F2, "not same # individuals for features %s and %s" % (keys[0], set)

    N = feature_sets[keys[0]].shape[0]
    inputs = np.zeros((N, 0))
    feature_names = []
    dim = {}
    dimsum = 0
    for set in keys:
        inputs_set = feature_sets[set].values
        dim[set] = inputs_set.shape[1]
        dimsum += dim[set]
        inputs = np.hstack((inputs, inputs_set))
        feature_names.extend(feature_sets[set].columns.tolist())

    if False:
        inputs.shape
        for j in keys: print(j + str(feature_sets[j].shape))
        import ipdb; ipdb.set_trace()

    #print "final size of inputs matrix is (%d, %d)" % inputs.shape
    return inputs, dim, dimsum, feature_names

#--------functions from featurization.py-----------
def featurize_data(data, learn_options, Y, gene_position, pam_audit=True, length_audit=True, quiet=True):
    '''
    assumes that data contains the 30mer
    returns set of features from which one can make a kernel for each one
    '''
    all_lens = data['30mer'].apply(len).values
    unique_lengths = np.unique(all_lens)
    num_lengths = len(unique_lengths)
    assert num_lengths == 1, "should only have sequences of a single length, but found %s: %s" % (num_lengths, str(unique_lengths))

    if not quiet:
        print("Constructing features...")
    t0 = time.time()

    feature_sets = {}

    if learn_options["nuc_features"]:
        # spectrum kernels (position-independent) and weighted degree kernels (position-dependent)
        get_all_order_nuc_features(data['30mer'], feature_sets, learn_options, learn_options["order"], max_index_to_use=30, quiet=quiet)

    check_feature_set(feature_sets)

    if learn_options["gc_features"]:
        gc_above_10, gc_below_10, gc_count = gc_features(data, length_audit)
        feature_sets['gc_above_10'] = pandas.DataFrame(gc_above_10)
        feature_sets['gc_below_10'] = pandas.DataFrame(gc_below_10)
        feature_sets['gc_count'] = pandas.DataFrame(gc_count)

    if learn_options["include_gene_position"]:
        # gene_position_columns = ["Amino Acid Cut position", "Percent Peptide", "Nucleotide cut position"]
        # gene_position_columns = ["Percent Peptide", "Nucleotide cut position"]

        for set in gene_position.columns:
            set_name = set
            feature_sets[set_name] = pandas.DataFrame(gene_position[set])
        feature_sets["Percent Peptide <50%"] = feature_sets["Percent Peptide"] < 50
        feature_sets["Percent Peptide <50%"]['Percent Peptide <50%'] = feature_sets["Percent Peptide <50%"].pop("Percent Peptide")

    if learn_options["include_gene_effect"]:
        print("including gene effect")
        gene_names = Y['Target gene']
        enc = sklearn.preprocessing.OneHotEncoder()
        label_encoder = sklearn.preprocessing.LabelEncoder()
        label_encoder.fit(gene_names)
        one_hot_genes = np.array(enc.fit_transform(label_encoder.transform(gene_names)[:, None]).todense())
        feature_sets["gene effect"] = pandas.DataFrame(one_hot_genes,
                                                       columns=["gene_%d" % i for i in range(one_hot_genes.shape[1])], index=gene_names.index)

    if learn_options['include_known_pairs']:
        feature_sets['known pairs'] = pandas.DataFrame(Y['test'])

    if learn_options["include_NGGX_interaction"]:
        feature_sets["NGGX"] = NGGX_interaction_feature(data, pam_audit)

    if learn_options["include_Tm"]:
        feature_sets["Tm"] = Tm_feature(data, pam_audit, learn_options=None)

    if learn_options["include_sgRNAscore"]:
        feature_sets["sgRNA Score"] = pandas.DataFrame(data["sgRNA Score"])

    if learn_options["include_drug"]:
        # feature_sets["drug"] = pandas.DataFrame(data["drug"])
        drug_names = Y.index.get_level_values('drug').tolist()
        enc = sklearn.preprocessing.OneHotEncoder()
        label_encoder = sklearn.preprocessing.LabelEncoder()
        label_encoder.fit(drug_names)
        one_hot_drugs = np.array(enc.fit_transform(label_encoder.transform(drug_names)[:, None]).todense())
        feature_sets["drug"] = pandas.DataFrame(one_hot_drugs, columns=["drug_%d" % i for i in range(one_hot_drugs.shape[1])], index=drug_names)

    if learn_options['include_strand']:
        feature_sets['Strand effect'] = (pandas.DataFrame(data['Strand']) == 'sense')*1

    if learn_options["include_gene_feature"]:
        feature_sets["gene features"] = gene_feature(Y, data, learn_options)

    if learn_options["include_gene_guide_feature"] > 0:
        tmp_feature_sets = gene_guide_feature(Y, data, learn_options)
        for key in tmp_feature_sets:
            feature_sets[key] = tmp_feature_sets[key]

    if learn_options["include_microhomology"]:
        feature_sets["microhomology"] = get_micro_homology_features(Y['Target gene'], learn_options, data)

    t1 = time.time()
    if not quiet:
        print("\t\tElapsed time for constructing features is %.2f seconds" % (t1-t0))

    check_feature_set(feature_sets)

    if learn_options['normalize_features']:
        assert("should not be here as doesn't make sense when we make one-off predictions, but could make sense for internal model comparisons when using regularized models")
        feature_sets = normalize_feature_sets(feature_sets)
        check_feature_set(feature_sets)

    return feature_sets


def check_feature_set(feature_sets):
    '''
    Ensure the # of people is the same in each feature set
    '''
    assert feature_sets != {}, "no feature sets present"

    N = None
    for ft in list(feature_sets.keys()):
        N2 = feature_sets[ft].shape[0]
        if N is None:
            N = N2
        else:
            assert N >= 1, "should be at least one individual"
            assert N == N2, "# of individuals do not match up across feature sets"

    for set in list(feature_sets.keys()):
        if np.any(np.isnan(feature_sets[set])):
            raise Exception("found Nan in set %s" % set)


def NGGX_interaction_feature(data, pam_audit=True):
    '''
    assuming 30-mer, grab the NGGX _ _ positions, and make a one-hot
    encoding of the NX nucleotides yielding 4x4=16 features
    '''
    sequence = data['30mer'].values
    feat_NX = pandas.DataFrame()
    # check that GG is where we think
    for seq in sequence:
        if pam_audit and seq[25:27] != "GG":
            raise Exception("expected GG but found %s" % seq[25:27])
        NX = seq[24]+seq[27]        
        NX_onehot = nucleotide_features(NX,order=2, feature_type='pos_dependent', max_index_to_use=2, prefix="NGGX")        
        # NX_onehot[:] = np.random.rand(NX_onehot.shape[0]) ##TESTING RANDOM FEATURE
        feat_NX = pandas.concat([feat_NX, NX_onehot], axis=1)
    return feat_NX.T


def get_all_order_nuc_features(data, feature_sets, learn_options, maxorder, max_index_to_use, prefix="", quiet=False):
    for order in range(1, maxorder+1):
        if not quiet:
            print("\t\tconstructing order %s features" % order)
        nuc_features_pd, nuc_features_pi = apply_nucleotide_features(data, order, learn_options["num_proc"],
                                                                     include_pos_independent=True, max_index_to_use=max_index_to_use, prefix=prefix)
        feature_sets['%s_nuc_pd_Order%i' % (prefix, order)] = nuc_features_pd
        if learn_options['include_pi_nuc_feat']:
            feature_sets['%s_nuc_pi_Order%i' % (prefix, order)] = nuc_features_pi
        check_feature_set(feature_sets)

        if not quiet:
            print("\t\t\t\t\t\t\tdone")


def countGC(s, length_audit=True):
    '''
    GC content for only the 20mer, as per the Doench paper/code
    '''
    if length_audit:
        assert len(s) == 30, "seems to assume 30mer"
    return len(s[4:24].replace('A', '').replace('T', ''))


def SeqUtilFeatures(data):
    '''
    assuming '30-mer'is a key
    get melting temperature features from:
        0-the 30-mer ("global Tm")
        1-the Tm (melting temperature) of the DNA:RNA hybrid from positions 16 - 20 of the sgRNA, i.e. the 5nts immediately proximal of the NGG PAM
        2-the Tm of the DNA:RNA hybrid from position 8 - 15 (i.e. 8 nt)
        3-the Tm of the DNA:RNA hybrid from position 3 - 7  (i.e. 5 nt)
    '''
    sequence = data['30mer'].values
    num_features = 1
    featarray = np.ones((sequence.shape[0], num_features))
    for i, seq in enumerate(sequence):
        assert len(seq) == 30, "seems to assume 30mer"
        featarray[i, 0] = SeqUtil.molecular_weight(str(seq))

    feat = pandas.DataFrame(pandas.DataFrame(featarray))
    return feat


def get_micro_homology_features(gene_names, learn_options, X):
    # originally was flipping the guide itself as necessary, but now flipping the gene instead

    print("building microhomology features")
    feat = pandas.DataFrame(index=X.index)
    feat["mh_score"] = ""
    feat["oof_score"] = ""

    #with open(r"tmp\V%s_gene_mismatches.csv" % learn_options["V"],'wb') as f:
    if True:
        # number of nulceotides to take to the left and right of the guide
        k_mer_length_left = 9
        k_mer_length_right = 21
        for gene in gene_names.unique():
            gene_seq = Seq.Seq(get_gene_sequence(gene)).reverse_complement()
            guide_inds = np.where(gene_names.values == gene)[0]
            print("getting microhomology for all %d guides in gene %s" % (len(guide_inds), gene))
            for j, ps in enumerate(guide_inds):
                guide_seq = Seq.Seq(X['30mer'][ps])
                strand = X['Strand'][ps]
                if strand=='sense':
                    gene_seq = gene_seq.reverse_complement()
                # figure out the sequence to the left and right of this guide, in the gene
                ind = gene_seq.find(guide_seq)
                if ind==-1:
                    gene_seq = gene_seq.reverse_complement()
                    ind = gene_seq.find(guide_seq)
                    #assert ind != -1, "still didn't work"
                    #print "shouldn't get here"
                else:
                    #print "all good"
                    pass
                #assert ind != -1, "could not find guide in gene"
                if ind==-1:
                    #print "***could not find guide %s for gene %s" % (str(guide_seq), str(gene))
                    #if.write(str(gene) + "," + str(guide_seq))
                    mh_score = 0
                    oof_score = 0
                else:
                    #print "worked"

                    assert gene_seq[ind:(ind+len(guide_seq))]==guide_seq, "match not right"
                    left_win = gene_seq[(ind - k_mer_length_left):ind]
                    right_win = gene_seq[(ind + len(guide_seq)):(ind + len(guide_seq) + k_mer_length_right)]

                    #if strand=='antisense':
                    #    # it's arbitrary which of sense and anti-sense we flip, we just want
                    #    # to keep them in the same relative alphabet/direction
                    #    left_win = left_win.reverse_complement()
                    #    right_win = right_win.reverse_complement()
                    assert len(left_win.tostring())==k_mer_length_left
                    assert len(right_win.tostring())==k_mer_length_right

                    sixtymer = str(left_win) + str(guide_seq) + str(right_win)
                    assert len(sixtymer)==60, "should be of length 60"
                    mh_score, oof_score = compute_score(sixtymer)

                feat.ix[ps,"mh_score"] = mh_score
                feat.ix[ps,"oof_score"] = oof_score
            print("computed microhomology of %s" % (str(gene)))

    return pandas.DataFrame(feat, dtype='float')


def local_gene_seq_features(gene_names, learn_options, X):

    print("building local gene sequence features")
    feat = pandas.DataFrame(index=X.index)
    feat["gene_left_win"] = ""
    feat["gene_right_win"] = ""

    # number of nulceotides to take to the left and right of the guide
    k_mer_length = learn_options['include_gene_guide_feature']
    for gene in gene_names.unique():
        gene_seq = Seq.Seq(get_gene_sequence(gene)).reverse_complement()
        for ps in np.where(gene_names.values==gene)[0]:
            guide_seq = Seq.Seq(X['30mer'][ps])
            strand = X['Strand'][ps]
            if strand=='sense':
                guide_seq = guide_seq.reverse_complement()
                #gene_seq = gene_seq.reverse_complement()
            # figure out the sequence to the left and right of this guide, in the gene
            ind = gene_seq.find(guide_seq)
            if ind ==-1:
                #gene_seq = gene_seq.reverse_complement()
                #ind = gene_seq.find(guide_seq)
                assert ind != -1, "could not find guide in gene"
            assert gene_seq[ind:(ind+len(guide_seq))]==guide_seq, "match not right"
            left_win = gene_seq[(ind - k_mer_length):ind]
            right_win = gene_seq[(ind + len(guide_seq)):(ind + len(guide_seq) + k_mer_length)]

            if strand=='antisense':
                # it's arbitrary which of sense and anti-sense we flip, we just want
                # to keep them in the same relative alphabet/direction
                left_win = left_win.reverse_complement()
                right_win = right_win.reverse_complement()
            assert not left_win.tostring()=="", "k_mer_context, %s, is too large" % k_mer_length
            assert not left_win.tostring()=="", "k_mer_context, %s, is too large" % k_mer_length
            assert len(left_win)==len(right_win), "k_mer_context, %s, is too large" % k_mer_length
            feat.ix[ps,"gene_left_win"] = left_win.tostring()
            feat.ix[ps,"gene_right_win"] = right_win.tostring()
        print("featurizing local context of %s" % (gene))

    feature_sets = {}
    get_all_order_nuc_features(feat["gene_left_win"], feature_sets, learn_options, learn_options["order"], max_index_to_use=sys.maxsize, prefix="gene_left_win")
    get_all_order_nuc_features(feat["gene_right_win"], feature_sets, learn_options, learn_options["order"], max_index_to_use=sys.maxsize, prefix="gene_right_win")
    return feature_sets

def gene_feature(Y, X, learn_options):
    '''
    Things like the sequence of the gene, the DNA Tm of the gene, etc.
    '''

    gene_names = Y['Target gene']

    gene_length = np.zeros((gene_names.values.shape[0], 1))
    gc_content = np.zeros((gene_names.shape[0], 1))
    temperature = np.zeros((gene_names.shape[0], 1))
    molecular_weight = np.zeros((gene_names.shape[0], 1))

    for gene in gene_names.unique():
        seq = get_gene_sequence(gene)
        gene_length[gene_names.values==gene] = len(seq)
        gc_content[gene_names.values==gene] = SeqUtil.GC(seq)
        temperature[gene_names.values==gene] = Tm.Tm_NN(seq, rna=False)
        molecular_weight[gene_names.values==gene] = SeqUtil.molecular_weight(seq, 'DNA')

    all = np.concatenate((gene_length, gc_content, temperature, molecular_weight), axis=1)
    df = pandas.DataFrame(data=all, index=gene_names.index, columns=['gene length',
                                                                     'gene GC content',
                                                                     'gene temperature',
                                                                     'gene molecular weight'])
    return df

def gene_guide_feature(Y, X, learn_options):
    #features, which are related to parts of the gene-local to the guide, and
    #possibly incorporating the guide or interactions with it

    #expensive, so pickle if necessary
    gene_file = r"..\data\gene_seq_feat_V%s_km%s.ord%s.pickle" % (learn_options['V'], learn_options['include_gene_guide_feature'], learn_options['order'])

    if False: #os.path.isfile(gene_file): #while debugging, comment out
        print("loading local gene seq feats from file %s" % gene_file)
        with open(gene_file, "rb") as f: feature_sets = pickle.load(f)
    else:
        feature_sets = local_gene_seq_features(Y['Target gene'], learn_options, X)
        print("writing local gene seq feats to file %s" % gene_file)
        with open(gene_file, "wb") as f: pickle.dump(feature_sets, f)

    return feature_sets


def gc_cont(seq):
    return (seq.count('G') + seq.count('C'))/float(len(seq))



def Tm_feature(data, pam_audit=True, learn_options=None):
    '''
    assuming '30-mer'is a key
    get melting temperature features from:
        0-the 30-mer ("global Tm")
        1-the Tm (melting temperature) of the DNA:RNA hybrid from positions 16 - 20 of the sgRNA, i.e. the 5nts immediately proximal of the NGG PAM
        2-the Tm of the DNA:RNA hybrid from position 8 - 15 (i.e. 8 nt)
        3-the Tm of the DNA:RNA hybrid from position 3 - 7  (i.e. 5 nt)
    '''

    if learn_options is None or 'Tm segments' not in list(learn_options.keys()):
        segments = [(19, 24), (11, 19), (6, 11)]
    else:
        segments = learn_options['Tm segments']

    sequence = data['30mer'].values
    featarray = np.ones((sequence.shape[0],4))

    for i, seq in enumerate(sequence):
        if pam_audit and seq[25:27]!="GG":
            raise Exception("expected GG but found %s" % seq[25:27])
        rna = False
        featarray[i,0] = Tm.Tm_NN(seq)        #30mer Tm
        featarray[i,1] = Tm.Tm_NN(seq[segments[0][0]:segments[0][1]]) #5nts immediately proximal of the NGG PAM
        featarray[i,2] = Tm.Tm_NN(seq[segments[1][0]:segments[1][1]])   #8-mer
        featarray[i,3] = Tm.Tm_NN(seq[segments[2][0]:segments[2][1]])      #5-mer

        #print "CRISPR"
        #for d in range(4):
        #    print featarray[i,d]
        #import ipdb; ipdb.set_trace()
    

    feat = pandas.DataFrame(featarray, index=data.index, columns=["Tm global_%s" % rna, "5mer_end_%s" %rna, "8mer_middle_%s" %rna, "5mer_start_%s" %rna])

    return feat

def gc_features(data, audit=True):
    gc_count = data['30mer'].apply(lambda seq: countGC(seq, audit))
    gc_count.name = 'GC count'
    gc_above_10 = (gc_count > 10)*1
    gc_above_10.name = 'GC > 10'
    gc_below_10 = (gc_count < 10)*1
    gc_below_10.name = 'GC < 10'
    return gc_above_10, gc_below_10, gc_count



def normalize_features(data,axis):
    '''
    input: Pandas.DataFrame of dtype=np.float64 array, of dimensions
    mean-center, and unit variance each feature
    '''
    data -= data.mean(axis)
    data /= data.std(axis)
    # remove rows with NaNs
    data = data.dropna(1)
    if np.any(np.isnan(data.values)): raise Exception("found NaN in normalized features")
    return data

def apply_nucleotide_features(seq_data_frame, order, num_proc, include_pos_independent, max_index_to_use, prefix=""):

    fast = True
    if include_pos_independent:
        feat_pd = seq_data_frame.apply(nucleotide_features, args=(order, max_index_to_use, prefix, 'pos_dependent'))
        feat_pi = seq_data_frame.apply(nucleotide_features, args=(order, max_index_to_use, prefix, 'pos_independent'))
        assert not np.any(np.isnan(feat_pd)), "nans here can arise from sequences of different lengths"
        assert not np.any(np.isnan(feat_pi)), "nans here can arise from sequences of different lengths"
        return feat_pd, feat_pi
    else:
        feat_pd = seq_data_frame.apply(nucleotide_features, args=(order, max_index_to_use, prefix, 'pos_dependent'))
        assert not np.any(np.isnan(feat_pd)), "found nan in feat_pd"
        return feat_pd

def get_alphabet(order, raw_alphabet = ['A', 'T', 'C', 'G']):
    alphabet = ["".join(i) for i in itertools.product(raw_alphabet, repeat=order)]
    return alphabet

def nucleotide_features(s, order, max_index_to_use, prefix="", feature_type='all', raw_alphabet = ['A', 'T', 'C', 'G']):
    '''
    compute position-specific order-mer features for the 4-letter alphabet
    (e.g. for a sequence of length 30, there are 30*4 single nucleotide features
          and (30-1)*4^2=464 double nucleotide features
    '''
    assert feature_type in ['all', 'pos_independent', 'pos_dependent']
    if max_index_to_use <= len(s):
        #print "WARNING: trimming max_index_to use down to length of string=%s" % len(s)
        max_index_to_use = len(s)

    if max_index_to_use is not None:
        s = s[:max_index_to_use]
    #assert(len(s)==30, "length not 30")
    #s = s[:30] #cut-off at thirty to clean up extra data that they accidentally left in, and were instructed to ignore in this way
    alphabet = get_alphabet(order, raw_alphabet = raw_alphabet)
    features_pos_dependent = np.zeros(len(alphabet)*(len(s)-(order-1)))
    features_pos_independent = np.zeros(np.power(len(raw_alphabet),order))

    index_dependent = []
    index_independent = []

    for position in range(0, len(s)-order+1, 1):
        for l in alphabet:
            index_dependent.append('%s%s_%d' % (prefix, l, position))

    for l in alphabet:
        index_independent.append('%s%s' % (prefix, l))

    for position in range(0, len(s)-order+1, 1):
        nucl = s[position:position+order]
        features_pos_dependent[alphabet.index(nucl) + (position*len(alphabet))] = 1.0
        features_pos_independent[alphabet.index(nucl)] += 1.0

        # this is to check that the labels in the pd df actually match the nucl and position
        assert index_dependent[alphabet.index(nucl) + (position*len(alphabet))] == '%s%s_%d' % (prefix, nucl, position)
        assert index_independent[alphabet.index(nucl)] == '%s%s' % (prefix, nucl)


    #index_independent = ['%s_pi.Order%d_P%d' % (prefix, order,i) for i in range(len(features_pos_independent))]
    #index_dependent = ['%s_pd.Order%d_P%d' % (prefix, order, i) for i in range(len(features_pos_dependent))]


    if np.any(np.isnan(features_pos_dependent)):
        raise Exception("found nan features in features_pos_dependent")
    if np.any(np.isnan(features_pos_independent)):
        raise Exception("found nan features in features_pos_independent")

    if feature_type == 'all' or feature_type == 'pos_independent':
        if feature_type == 'all':
            res = pandas.Series(features_pos_dependent,index=index_dependent), pandas.Series(features_pos_independent,index=index_independent)
            assert not np.any(np.isnan(res.values))
            return res
        else:
            res = pandas.Series(features_pos_independent, index=index_independent)
            assert not np.any(np.isnan(res.values))
            return res

    res = pandas.Series(features_pos_dependent, index=index_dependent)
    assert not np.any(np.isnan(res.values))    
    return res

def nucleotide_features_dictionary(prefix=''):
    seqname = ['-4', '-3', '-2', '-1']
    seqname.extend([str(i) for i in range(1,21)])
    seqname.extend(['N', 'G', 'G', '+1', '+2', '+3'])

    orders = [1, 2, 3]
    sequence = 30
    feature_names_dep = []
    feature_names_indep = []
    index_dependent = []
    index_independent = []

    for order in orders:
        raw_alphabet = ['A', 'T', 'C', 'G']
        alphabet = ["".join(i) for i in itertools.product(raw_alphabet, repeat=order)]
        features_pos_dependent = np.zeros(len(alphabet)*(sequence-(order-1)))
        features_pos_independent = np.zeros(np.power(len(raw_alphabet),order))

        index_dependent.extend(['%s_pd.Order%d_P%d' % (prefix, order, i) for i in range(len(features_pos_dependent))])
        index_independent.extend(['%s_pi.Order%d_P%d' % (prefix, order,i) for i in range(len(features_pos_independent))])

        for pos in range(sequence-(order-1)):
            for letter in alphabet:
                feature_names_dep.append('%s_%s' % (letter, seqname[pos]))

        for letter in alphabet:
            feature_names_indep.append('%s' % letter)

        assert len(feature_names_indep) == len(index_independent)
        assert len(feature_names_dep) == len(index_dependent)

    index_all = index_dependent + index_independent
    feature_all = feature_names_dep + feature_names_indep

    return dict(list(zip(index_all, feature_all)))

def normalize_feature_sets(feature_sets):
    '''
    zero-mean, unit-variance each feature within each set
    '''

    print("Normalizing features...")
    t1 = time.time()

    new_feature_sets = {}
    for set in feature_sets:
         new_feature_sets[set] = normalize_features(feature_sets[set],axis=0)
         if np.any(np.isnan(new_feature_sets[set].values)):
             raise Exception("found Nan feature values in set=%s" % set)
         assert new_feature_sets[set].shape[1] > 0, "0 columns of features"
    t2 = time.time()
    print("\t\tElapsed time for normalizing features is %.2f seconds" % (t2-t1))

    return new_feature_sets

#-------util.py functions------

def get_gene_sequence(gene_name):
    try:
        gene_file = '../../gene_sequences/%s_sequence.txt' % gene_name
        #gene_file = '../gene_sequences/%s_sequence.txt' % gene_name
        #gene_file = 'gene_sequences/%s_sequence.txt' % gene_name
        with open(gene_file, 'rb') as f:
            seq = f.read()
            seq = seq.replace('\r\n', '')
    except:
        raise Exception("could not find gene sequence file %s, please see examples and generate one for your gene as needed, with this filename" % gene_file)

    return seq

#-----microhomology.py-----
def compute_score(seq, tmpfile1="1.before removing duplication.txt", tmpfile2="2.all microhomology patterns.txt", verbose=False):
    length_weight=20.0 
    left=30        # Insert the position expected to be broken. 
    right=len(seq)-int(left) 
    #print 'length of seq = '+str(len(seq)) 
     
    file_temp=open(tmpfile1, "w") 
    for k in range(2,left)[::-1]: 
            for j in range(left,left+right-k+1): 
                    for i in range(0,left-k+1): 
                            if seq[i:i+k]==seq[j:j+k]: 
                                    length=j-i 
                                    file_temp.write(seq[i:i+k]+'\t'+str(i)+'\t'+str(i+k)+'\t'+str(j)+'\t'+str(j+k)+'\t'+str(length)+'\n') 
    file_temp.close() 
     
    ### After searching out all microhomology patterns, duplication should be removed!! 
    f1=open(tmpfile1, "r") 
    s1=f1.read() 
     
    f2=open(tmpfile2, "w") #After removing duplication 
    f2.write(seq+'\t'+'microhomology\t'+'deletion length\t'+'score of a pattern\n') 
     
    if s1!="": 
            list_f1=s1.strip().split('\n') 
            sum_score_3=0 
            sum_score_not_3=0 
     
            for i in range(len(list_f1)): 
                    n=0 
                    score_3=0 
                    score_not_3=0 
                    line=list_f1[i].split('\t') 
                    scrap=line[0] 
                    left_start=int(line[1]) 
                    left_end=int(line[2]) 
                    right_start=int(line[3]) 
                    right_end=int(line[4]) 
                    length=int(line[5]) 
     
                    for j in range(i): 
                            line_ref=list_f1[j].split('\t') 
                            left_start_ref=int(line_ref[1]) 
                            left_end_ref=int(line_ref[2]) 
                            right_start_ref=int(line_ref[3]) 
                            right_end_ref=int(line_ref[4]) 
     
                            if (left_start >= left_start_ref) and (left_end <= left_end_ref) and (right_start >= right_start_ref) and (right_end <= right_end_ref): 
                                    if (left_start - left_start_ref)==(right_start - right_start_ref) and (left_end - left_end_ref)==(right_end - right_end_ref): 
                                            n+=1 
                            else: pass 
                           
                    if n == 0: 
                            if (length % 3)==0: 
                                    length_factor = round(1/exp((length)/(length_weight)),3) 
                                    num_GC=len(findall('G',scrap))+len(findall('C',scrap)) 
                                    score_3=100*length_factor*((len(scrap)-num_GC)+(num_GC*2)) 
                                     
                            elif (length % 3)!=0: 
                                    length_factor = round(1/exp((length)/(length_weight)),3) 
                                    num_GC=len(findall('G',scrap))+len(findall('C',scrap)) 
                                    score_not_3=100*length_factor*((len(scrap)-num_GC)+(num_GC*2)) 
     
                            f2.write(seq[0:left_end]+'-'*length+seq[right_end:]+'\t'+scrap+'\t'+str(length)+'\t'+str(100*length_factor*((len(scrap)-num_GC)+(num_GC*2)))+'\n') 
                    sum_score_3+=score_3 
                    sum_score_not_3+=score_not_3 
     
            mh_score = sum_score_3+sum_score_not_3
            oof_score = (sum_score_not_3)*100/(sum_score_3+sum_score_not_3)
            if verbose:
                print('Microhomology score = ' + str(mh_score)) 
                print('Out-of-frame score = ' + str(oof_score)) 
    f1.close() 
    f2.close()
    return mh_score, oof_score
🚀 Quickstart 
==============

Installation
**************
PEGG is available through the python package index. To install, use pip: 

.. code-block:: python

   pip install pegg

Required Input Files
****************************
In order to use PEGG, you need at minimum 2 input files:

(1) A dataframe containing a list of input mutations

(2) A reference genome to work from.

An optional, but highly reccomended input file for use with the reference genome build GrCh37 is:

(3) chrom_dict, which provides on/off target scores for each protospacer occuring in exonic regions of GrCh37.

Reference files for (2) and (3), as well as a sample mutation dataset (1), are available at the following link:

- `Reference Files Dropbox Link <https://www.dropbox.com/sh/h6fdvpv3tyny27q/AADYVOkJe12XZiD4pf3_WXuga?dl=0>`_

Reference File Formatting
**************************
PEGG requires that the reference files are in a specific format. This should be the only painful part of setting up your environment
to allow PEGG to generate pegRNAs. 

(1) Input mutations
~~~~~~~~~~~~~~~~~~~~~

First let's load in a dataframe containing a list of input mutations (1):

.. code-block:: python

   import pegg
   import pandas as pd

   #------(1) loading in input mutatations-------------
   filepath = '.../2020-06-16-MSK-IMPACT_EDITED.txt'
   mutant_input = pd.read_csv(filepath, sep='\t')

The above codeblock shows how to load in the sample mutation input file available in the `dropbox of reference files <https://www.dropbox.com/sh/h6fdvpv3tyny27q/AADYVOkJe12XZiD4pf3_WXuga?dl=0>`_.
To use **your own mutation input**, the input mutations must follow the format shown in the dataframe below:

.. image:: mutant_input.png

Namely, 7 columns are required, providing information about the input mutations relative to the provided reference genome.
The input mutation dataframe **must match the column naming scheme shown above**:

1. Hugo_Symbol: the name of the gene within which the mutations falls. For non-coding locations, include some identifer.

2. Chromosome: which chromosome the mutation falls in. Use integers or 'X' and 'Y'.

3. Start_Position: start position of mutation in the reference genome (mutations should only be reported on the + strand)

4. End_Position: end position of mutation in the reference genome

5. Variant_Type: what type of mutation is it. Options: "SNP", "ONP", "INS", "DEL"

6. Reference_Allele: what is the reference allele. For insertions, this can be set to "-".

7. Tumor_Seq_Allele2: what is the mutant allele (i.e. what is the mutation sequence). For deletions, this can be set to "-".

(1a) Compatible Input Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All of the mutation files listed in the cBioPortal are in the correct format for input into PEGG.
They can be accessed at the following link: http://www.cbioportal.org/datasets

See the jupyter notebook tutorial for a more in depth case study of using these reference databases.


(1b) ClinVar variants
~~~~~~~~~~~~~~~~~~~~~~~
PEGG also includes a built in tool that allows users to put a list of ClinVar variants in the correct format for interpretation by PEGG. 
This is done using the **pegg.clinvar_VCF_translator()** function, where users provide a clinvar.vcf.gz file and a list of Variation IDs that 
correspond with the variants they want to translate. These variation IDs are the identifiers for ClinVar variants:

.. image:: var_ids.png

See the below codeblock for the precise syntax:

.. code-block:: python

   #this is the filepath to the .vcf.gz file (needs to be updated according to user)
   filepath = '.../clinvar_20221001.vcf.gz'
   variation_ids = [925574, 925434, 926695, 925707, 325626, 1191613,308061] #list of ClinVar variants to translate
   clinvar = pegg.clinvar_VCF_translator(filepath, variation_ids)

This outputs a dataframe in the correct format for generating pegRNAs using PEGG:

.. image:: out_clinvar.png

ClinVar VCF files can be accessed here: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/

A sample vcf.gz file with GrCh37 coordinates is provided in the `dropbox of reference files <https://www.dropbox.com/sh/h6fdvpv3tyny27q/AADYVOkJe12XZiD4pf3_WXuga?dl=0>`_ as well.

(2) Reference Genome
~~~~~~~~~~~~~~~~~~~~~

And now let's load in the reference genome. I suggest using the codeblock below if you're using GrCh37 or GrCh38 as a reference file.
Otherwise, the format needs to be as follows:

- **records** = a list containing the chromosome sequences (preferably parsed by SeqIO).

- **index_list** = a list containing the indeces in **records** corresponding to chromosome 1-22, and chrX and chrY.

   - For example, **records[index_list[0]] should refer to chromosome 1** (numbering starts at 0, not 1 in index_list due to python conventions)
   - And **records[index_list[1]] should refer to chromosome 2**, etc.
   - For the X and Y chromosome, records[index_list[22]] should refer to chrX, and records[index_list[23]] should refer to chrY

In other words, **records** is a list containing the chromosome sequences, but likely also other unwanted sequences (alternate/unlocalized sequences), 
while **index_list** provides the indeces needed to place the chromosome sequences in records in order from chromosome 1-22, and chrX and chrY.

- An alternate way to format these input files is to have **records** list have the chromosome sequences in order from 1-22, chrX, and chrY, and set **index_list** = list(range(24))

See the following codeblock for how to set up **records** and **index_list** using the `provided reference genome files, GrCh37 & GrCh38 <https://www.dropbox.com/sh/h6fdvpv3tyny27q/AADYVOkJe12XZiD4pf3_WXuga?dl=0>`_.

.. code-block:: python

   #filepath to .gz containing reference genome (here = GrCh37)
   filepath = './GCF_000001405.25_GRCh37.p13_genomic.fna.gz'

   records, index_list = pegg.genome_loader(filepath)

(3) On/Off Target scores
~~~~~~~~~~~~~~~~~~~~~~~~~
The last reference file (3) is chrom_dict, which provides on/off target scores for protospacers in exonic regions of GrCh37.
These scores are used in the ranking and filtration of pegRNAs. Curently this is only available for GrCh37, but a chrom_dict reference file for
GrCh38 will be provided shortly and available in the `dropbox of reference files <https://www.dropbox.com/sh/h6fdvpv3tyny27q/AADYVOkJe12XZiD4pf3_WXuga?dl=0>`_.

To import this reference file, use the following syntax:

.. code-block:: python

   #---------loading in on/off-target efficiencies of sgRNAs in exonic regions of GrCh37---------------------
   file = '.../chrom_dict.pickle'
   chrom_dict = pd.read_pickle(file)

If you are using an alternative genome, or don't wish to use chrom_dict, set it equal to 'none':


.. code-block:: python

   chrom_dict = 'none'


Using PEGG
***********
Now that our reference files are loaded in, and PEGG is imported as a module, using PEGG is simple.
We simply need to specify parameters which correspond with the different options associated with prime editing, 
as depicted in the visualization below:

.. image:: PE_schematic.png


Namely, the user must specify:

1. Select mutations within mutant_input to generate pegRNAs for. This is done by providing a list of indeces that correspond with the desired mutations. The alternative is simply setting this to **all mutations in the datasets, by having mut_idx_list = list(range(len(mutant_input)))**.

2. PAM sequence (string format)

3. How many guides to return per mutation

4. A list of RTT and PBS lengths to search.

.. code-block:: python
   
   #specify which mutations within mutant_input to generate pegRNAs for
   #here we're going to just generate pegRNAs for 1 mutation, corresponding to row 4 of mutant_input
   mut_idx_list = [4] 
   PAM = 'NGG' 
   guides_per_mut = 5  #specify how many pegRNAs to return for each mutation
   RTT_lengths = [20,25,30] #specify RTT lengths and PBS lengths to search
   PBS_lengths = [5,7,10]
   minus_seqs = pegg.minus_seq_generator(records, index_list)

   #now generating the pegRNAs
   run_output = pegg.run(mutant_input, mut_idx_list, records, index_list, minus_seqs, chrom_dict, PAM, RTT_lengths, PBS_lengths, guides_per_mut)
   

Visualization Tools
********************

PEGG has built in tools for visualizing the pegRNAs it generates, providing the ability to spot-check designs.

In the sample below, we generate our pegRNAs using the run() function and then select a pegRNA from the resulting
output dataframe to visualize, using **pegg.pegrna_display()**:


.. code-block:: python

   pegRNA_df_loc=0 #choosing which guide to display from the dataframe
   h = pegg.pegrna_display(run_output, pegRNA_df_loc, records, index_list)

.. image:: pegviz.png

There's another built-in tool for visualizing the 3' extension (RTT and PBS sequence) of pegRNAs.
In the example below, we use it to visualize the 3' extensions of the first 4 guides in the output using
**pegg.align_display()**:

.. code-block:: python

   pegg.align_display(run_output[0:4], records, index_list)

.. image:: align_display.png

Oligo Generation
*****************

To automatically generate oligonucleotides that contain the pegRNAs designed using PEGG, the **pegg.oligo_generator()**
function provides multiple options, and produces both a **pegRNA oligo** and an **epegRNA oligo** (with a 3' structural motif, `tevopreQ1 <https://www.nature.com/articles/s41587-021-01039-7>`_).


A unique feature of PEGG is the option to include a `sensor region <https://www.nature.com/articles/s41587-021-01172-3>`_  in the oligo. 
This sensor region is a synthetic version of the endogenous target site, providing the ability to measure a proxy of editing outcomes at the endogenous locus.
This approach can be useful in the context of a library of pegRNAs, allowing for the measurement of pegRNA enrichment/depletion *as well as* a proxy of editing outcomes
with a single NGS amplicon. The below schematic shows a schematic of the oligos that are output with sensor=True or sensor=False:

.. image:: oligos.png

Additionally, users need to specify whether they want to append a 'G' nucleotide to the beginning of the protospacer. 
This is reccomended in the original Anzalone et al., 2019 prime editing paper. The sensor and append_proto_G options are both set to True in the below example.

.. code-block:: python

   oligos_w_sensor = pegg.oligo_generator(run_output, append_proto_G=True, sensor=True)


This returns a dataframe that has the oligos appended as columns ('pegRNA_oligo' and 'epegRNA_tevopreQ1_oligo' are the column names).

Users can also specify which 3' and 5' adapter sequences they want to use, or simply leave these options blank
and use the built-in adapters provided by the authors. In addition, users can specify to use a different gRNA scaffold,
or use the canonical gRNA scaffold provided by the authors. In the above example, these parameters 
(3_prime_adapter, 5_prime_adapter, and gRNA_scaff) are left empty, so the default versions provided by the author are used.

See the complete function documentation tab for more information about **pegg.oligo_generator()**.


Library Generation
********************
PEGG also includes automated library generation and visualization functions.
These provide the ability to automatically select all of the mutations associated with a particular gene, 
generate pegRNAs for these mutations, and generate neutral pegRNAs that introduce silent mutations as internal controls.

The code below shows how to generate the neutral/silent substitutions based on inputting information about a gene
as well as providing a list of the coding sequence locations of the relevant transcript. This list is generated manually in the example 
below. The jupyter notebook tutorial shows how this step can be automated based on using available gene annotations.

.. code-block:: python

   gene_name='TP53'
   strand = '-'
   chrom='chr17'
   #listing CDS of transcript ordered by +end
   start_end_cds = [[7572930, 7573008],
   [7573927, 7574033],
   [7576853, 7576926],
   [7577019, 7577155],
   [7577499, 7577608],
   [7578177, 7578289],
   [7578371, 7578554],
   [7579312, 7579590],
   [7579700, 7579721],
   [7579839, 7579912]]
   neutral_p53 = pegg.neutral_substitutions(gene_name, chrom, strand, start_end_cds, records, index_list)

This generates a dataframe of all possible neutral mutations:

.. image:: neutral_sub.png

The above function is actually not needed to generate these libraries with internal controls included.
This can be done by simply running the below function: 

.. code-block:: python

   control_fraction=.01 #what fraction of the library do you want to be neutral/silent pegRNAs
   library_input = pegg.library_input_generator(mutant_input, gene_name, chrom, strand, start_end_cds, records, index_list, control_fraction)

Once this library input is generated, this can simply be fed into the **pegg.run()** function as shown previously.
In addition, there are built in library visualization tools. To use these, the user needs to add some information back into
the library_input dataframe. Namely, neutral guides need to be labelled, and HGVSp information needs to be added back to the dataframe
if it's available:

.. code-block:: python

   #generating the pegRNA library
   #same input required as shown previously
   ranked_filtered = pegg.run(mutant_input, mut_idx_list, records, index_list, minus_seqs, chrom_dict, PAM, RTT_lengths, PBS_lengths, guides_per_mut)

   #adding HGVSp information back to the dataframe if it's available...
   hg = []
   for i, val in ranked_filtered.iterrows():
      idx = val['mutant index']
      hgvsp = mutant_input.loc[[idx]]['HGVSp'].values[0]
      hg.append(hgvsp)
      
   #also add in information for identifying neutral mutations
   class_mut = []
   for i, val in ranked_filtered.iterrows():
      idx = val['mutant index']
      neut = mutant_input.loc[[idx]]['classification'].values[0]
      class_mut.append(neut)

   ranked_filtered['HGVSp']=hg
   ranked_filtered['classification']=class_mut

Once this is done, the libraries can be visualized using either of the two functions below:

.. code-block:: python

   pegg.lollipop_library(ranked_filtered, gene_name, start_end_cds, strand, plot=True)


.. image:: lollipop.png


.. code-block:: python

   pegg.matrix_rep_library(ranked_filtered, gene_name, start_end_cds, strand, plot=True)

.. image:: matrix_lib.png

More information about the library generation functionality is provided in the jupyter notebook tutorial.

Jupyter Notebook Tutorial
**************************
A jupyter notebook version of the PEGG tutorial can be accessed at the following link: 

`Jupyter Notebook Tutorial <https://github.com/samgould2/PEGG/blob/main/examples/PEGG_example.ipynb>`_


A Note on RAM
**************
Importing a reference genome into the local environment requires ~4 Gb of RAM. Chrom_dict is also a large file.
It's reccomended to use a machine with *at least*  16 Gb of RAM, though more is preferable, when running pegg.


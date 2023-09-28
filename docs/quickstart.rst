🚀 Quickstart 
==============

Installation
**************
PEGG is available through the python package index. To install, use pip: 

.. code-block:: python

   pip install pegg

To modify and/or download specific python files, download the package from the `Github Repository <https://github.com/samgould2/PEGG2.0>`_ .

In order to use PEGG to design pegRNAs, you must provide a set of mutations in a pandas DataFrame for PEGG to use.
There are 3 acceptable input formats that PEGG accepts

(1) WT-ALT Format
*******************
In this format, the user feeds in a set of wildtype sequences, and the corresponding set of desired mutant sequences.
PEGG then uses a built-in aligner to determine the mutation from these sequences, and design pegRNAs for these variants.
Note: Complex INDELs will break this function because they are difficult to align. If you want to design INDELs, use one of the two formats below.

The pandas DataFrame must have at least 2 columns, with the following column names **(1) 'WT'** and **(2) 'ALT'**

(2) PrimeDesign Format
************************

In this format, users feed in a set of sequences with the desired edit in `PrimeDesign <https://primedesign.pinellolab.partners.org/>`_ format.
For example, an A>G SNP would look like: "AATCG(A/G)GCTAG", a "GTT" insertion would be: "AATCG(/GTT)GCTAG", a "GTT" deletion: "AATCG(GTT/)GCTAG", and a "C>GTT" INDEL: "AATCG(C/GTT)GCTAG".
It is reccomended that you include at least 100 nt of context sequence on either side of the variant.

The pandas DataFrame must have at least 1 column, with the following column name: **(1) 'SEQ'** .
This column contains the PrimeDesign formatted sequences.


(3) cBioPortal Format
***********************

In this format, users simply provide a set of variants with their genome coordinates and reference and alternate alleles, as well as the variant type.
This format is compatible with all of the datasets in the `cBioPortal <http://www.cbioportal.org/datasets>`_ , which contains cancer-associated variant datasets.

If you want to build you own mutation dataset, it must have the following columns, with the associated column header names, in .csv format:

1. Chromosome: which chromosome the mutation falls in. Use integers or 'X' and 'Y'.

2. Start_Position: start position of mutation in the reference genome (mutations should only be reported on the + strand)

3. End_Position: end position of mutation in the reference genome

4. Variant_Type: what type of mutation is it. Options: "SNP", "ONP", "INS", "DEL"

5. Reference_Allele: what is the reference allele. For insertions, this can be set to "-".

6. Tumor_Seq_Allele2: what is the mutant allele (i.e. what is the mutation sequence). For deletions, this can be set to "-".

An example for loading in a dataset is provided here:

.. code-block:: python

   import pandas as pd

   #------(1) loading in input mutatations-------------
   filepath = '.../2020-06-16-MSK-IMPACT_EDITED.txt'
   mutant_input = pd.read_csv(filepath, sep='\t')

.. image:: mutant_input.png

(3b) Reference Genome
***********************

In addition, to use the "cBioPortal" input format, you must provide a genome that PEGG can use to find the associated sequences. 
There is a built in genome loader function that you can use, or you can format the chromosome sequences in a dictionary, with the keys as the chromosome identifiers:

.. code-block:: python

   from pegg.prime import pegg2 

   #filepath to .gz containing reference genome (here = GrCh37)
   filepath = './GCF_000001405.25_GRCh37.p13_genomic.fna.gz'

   chrom_dict, i = pegg2.genome_loader(filepath)

You can access the genome files at this link: `Reference Files Dropbox Link <https://www.dropbox.com/sh/5xsdzyiyrjiu9pf/AADiFFA3BQ3vX7swja-i2NBqa?dl=0>`_

It has been tested on human and mouse genomes.

(4) ClinVar mutations
***********************

PEGG also has a built-in function for translating ClinVar IDs to the "cBioPortal" input format. To do so, download a ClinVar vcf.gz file,
and choose your desired Variation IDs that you wish to model. These vcf.gz files are available in the `Reference Files Dropbox Link <https://www.dropbox.com/sh/5xsdzyiyrjiu9pf/AADiFFA3BQ3vX7swja-i2NBqa?dl=0>`_ (under the clinvar folder)
or the more up to date files can be accessed here: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/

These variation IDs are the identifiers for ClinVar variants:

.. image:: var_ids.png

See the below codeblock for the precise syntax:

.. code-block:: python

   from pegg.prime import pegg2 

   #filepath to the clinvar vcf.gz file
   filepath = '.../GrCh37_clinvar_20230923.vcf.gz'
   variation_ids = [925574, 925434, 926695, 925707, 325626, 1191613, 308061, 361149, 1205375, 208043]
   df = pegg2.clinvar_VCF_translator(filepath, variation_ids)


Generating pegRNAs
********************



pegRNA Design Options
************************



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


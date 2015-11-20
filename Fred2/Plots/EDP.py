from Fred2.Core import Protein, Generator, Allele, Peptide
import Fred2.Core.Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
import Fred2.Core.Allele
import Fred2.EpitopePrediction.PSSM as PSSM_df
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import re


def epitope_density_plot(epitopeResult, transcript_id, method=None):
    """
        Description
    """
    # -----------------------------------------------------------------------------------
    # what protein is the peptide from
    # -----------------------------------------------------------------------------------

    s = epitopeResult
    print(s)
    peptide_index = epitopeResult.index.levels[0]
    protein_seq = peptide_index[0].get_protein(transcript_id)

    # -----------------------------------------------------------------------------------
    # Create peptide_per_position dictionary
    # -----------------------------------------------------------------------------------

    peptide_per_pos = {}
    for i in xrange(0,len(peptide_index)):
        f = peptide_index[i].get_protein_positions(transcript_id)
        peptide_per_pos[str(peptide_index[i])] = f
    print(peptide_per_pos)

    # -----------------------------------------------------------------------------------
    # Initialize and fill epitope_matrix with values
    # -----------------------------------------------------------------------------------

    epitope_matrix = np.zeros(shape=(len(epitopeResult.columns), len(protein_seq)))
    print(epitope_matrix)

    # -----------------------------------------------------------------------------------
    # define plot
    # -----------------------------------------------------------------------------------
    #plt.rc('font', family='','lines ,lw=2 ,c='')


    cmap = plt.cm.get_cmap('Reds')
    cmap.set_bad(color='0.75', alpha=None)

    fig = plt.figure()
    gs = matplotlib.gridspec.GridSpec(1, 1)
    ax1 = plt.subplot(gs[0])
    #ax1 = fig.add_subplot(gs[0])
    ax1.plot([1, 2])
    plt.title('Epitope_Density_Plot')
    #plt.show()

# -----------------------------------------------------------------------------------
# Create test case
# -----------------------------------------------------------------------------------

protein = Protein("ASDERWQTGHKILPMNVFCY", gene_id=1, transcript_id="someID")
peps = Generator.generate_peptides_from_proteins(protein, 9)
result = EpitopePredictorFactory("BIMAS").predict(peps, alleles=Allele("HLA-A*02:01"))

epitope_density_plot(result, "someID", method="BIMAS")





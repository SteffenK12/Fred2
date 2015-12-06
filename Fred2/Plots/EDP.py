from Fred2.Core import Protein, Generator, Allele, Peptide
import Fred2.Core.Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from pylab import rcParams
import Fred2.Core.Allele
import Fred2.EpitopePrediction.PSSM as PSSM_df
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import re


def epitope_density_plot(epitopeResult, transcript_id, method=None, cm_max=None, cm_min=None):
    """
        Description

        :param epitopeResult: Pandas Dataframe of an epitope prediction as :class:`~Fred2.Core.Result.EpitopePredictionResult` object
        :type epitopeResult: :class:`~Fred2.Core.Base.AEpitopePredictionResult`
        :param str transcript_id: The unique transcript ID of the :class:`~Fred2.Core.Protein.Protein` used in the plot
        :param method: Method used for epitope prediction from :class:`~Fred2.EpitopePrediction.PSSM.APSSMEpitopePrediction` resulting in an  :class:`~Fred2.Core.Result.EpitopePredictionResult` object
        :type method:
        :param int cm_max: Maximum density for epitope_density_plot, default=None
        :param int cm_min: Minimum density for epitope_density_plot, default=None
        :results: Draws an epitope_density_plot for a given EpitopePredictionResult object and saves it as epitope_density_plot.png
        :rtype: Figure object of the created epitope_density_plot
    """
    # -----------------------------------------------------------------------------------
    # what protein is the peptide from
    # -----------------------------------------------------------------------------------

    peptide_index = epitopeResult.index.levels[0]
    protein_seq =None

    # -----------------------------------------------------------------------------------
    # Create peptide_per_position dictionary
    # -----------------------------------------------------------------------------------

    peptide_per_pos = {}
    for i in peptide_index:
        if transcript_id in i.proteins:
            if protein_seq is None:
                protein_seq = i.proteins[transcript_id]
            peptide_per_pos[i] = i.get_protein_positions(transcript_id)


    # -----------------------------------------------------------------------------------
    # Initialize and fill epitope_matrix with values
    # -----------------------------------------------------------------------------------

    epitope_matrix = np.zeros((len(epitopeResult.columns), len(protein_seq)))
    #for i in xrange(0, len(peptide_index)):
    #    key = peptide_index[i]
    #    print(key)
    for a, allele in enumerate(epitopeResult.columns):
        for key, values in peptide_per_pos.iteritems():
            x = peptide_per_pos.get(key)
            for i in values:
                for k in xrange(i,i+8):
                    epitope_matrix[a][k] += epitopeResult.loc[(key, epitopeResult.index.levels[1][0]), allele]
    print(epitope_matrix)

    # -----------------------------------------------------------------------------------
    # define plot
    # -----------------------------------------------------------------------------------
    #plt.rc('font', family='','lines ,lw=2 ,c='')


    cmap = plt.cm.get_cmap('Reds')
    cmap.set_bad(color='0.75', alpha=None)


    fig = plt.figure()
    #fig.set_size_inches(epitope_matrix.shape[1]*0.5,epitope_matrix.shape[0], forward=True)
    gs = matplotlib.gridspec.GridSpec(1, 1, height_ratios=[1], width_ratios=[1])
    ax1 = plt.subplot(gs[0])
    if cm_max and cm_min is None:
        cm = plt.pcolormesh(epitope_matrix, cmap=cmap, vmax=epitope_matrix.max(), vmin=epitope_matrix.min())
    else:
        cm = plt.pcolormesh(epitope_matrix, cmap=cmap, vmax=cm_max, vmin=cm_min)

    #ax1 = fig.add_subplot(gs[0])
    #plt.xlim([0, len(protein_seq)])

    LABEL_SIZE = 10
    LABEL_X_OFFSET = 0.5
    LABEL_Y_OFFSET = 0.5
    ax1.xaxis.tick_top()
    ax1.set_xticks(np.arange(epitope_matrix.shape[1]))
    ax1.set_yticks(np.arange(epitope_matrix.shape[0]))
    y_label_distance = 0.7
    x_label_distance = 0.375
    column_labels = list(protein_seq)
    for i in range(epitope_matrix.shape[0]):
        # write yaxis labels as text
        ax1.text(-2.25,y_label_distance, str(epitopeResult.columns.values[i]), size = 10)
        y_label_distance += 1
    for i in range(epitope_matrix.shape[1]):
        # write xaxis labels as text
        ax1.text(x_label_distance, -0.275, str(column_labels[i]), size = 10)
        x_label_distance += 1
    plt.tick_params(axis='x', labelsize=LABEL_SIZE)
    plt.tick_params(axis='y', labelsize=LABEL_SIZE)
    ax1.set_yticklabels('')
    ax1.set_xticklabels('')
    ax1.set_ylim(20, 0)
    #fig.savefig('epitope_density_plot.pdf') #save figure with given filename
    plt.grid(b=True, which='major', axis='both', color='black', linestyle='-')
    cbar = fig.colorbar(cm)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('epitope probability', rotation=270)
    plt.show()
    plt.close(fig)

# -----------------------------------------------------------------------------------
# Create test case
# -----------------------------------------------------------------------------------

protein = Protein("ASDERWQTGHKILPMNVFCY", gene_id=1, transcript_id="someID")
peps = Generator.generate_peptides_from_proteins(protein, 9)
result = EpitopePredictorFactory("BIMAS").predict(peps, alleles=[Allele("HLA-A*02:01"), Allele("C*07:02"), Allele("B*27:02"), Allele("B*39:01"), Allele("B*51:03"), Allele("B*40:06"), Allele("B*38:01")])

epitope_density_plot(result, "someID", method="BIMAS")





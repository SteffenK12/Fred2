from __future__ import division
from Fred2.EpitopePrediction import EpitopePredictorFactory
from pylab import rcParams
import Fred2.Core.Allele
import Fred2.EpitopePrediction.PSSM as PSSM_df
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.gridspec as gridspec
import math
import operator



def epitope_cluster_plot(epitopeResult, transcript_id, method=None, cm_max=None, cm_min=None, cmp=lambda a,b: True, threshold=None,  savepath='output.pdf'):
    """

        :param epitopeResult: Pandas Dataframe of an epitope prediction
        as :class:`~Fred2.Core.Result.EpitopePredictionResult` object
        :type epitopeResult: :class:`~Fred2.Core.Base.AEpitopePredictionResult`
        :param str transcript_id: The unique transcript ID of the :class:`~Fred2.Core.Protein.Protein` used in the plot
        :param method: Method used for epitope prediction
        from :class:`~Fred2.EpitopePrediction.PSSM.APSSMEpitopePrediction` resulting in
        an :class:`~Fred2.Core.Result.EpitopePredictionResult` object
        :type method:
        :param int cm_max: Maximum density for epitope_cluster_plot, default=None
        :param int cm_min: Minimum density for epitope_cluster_plot, default=None
        :param cmp: Comparator used to check predicted value vs threshold
        :param threshold: Threshold used to filter epitopes
        :param savepath: Path and name of the file to be saved
        :results: Draws an epitope_cluster_plot for a given EpitopePredictionResult object and saves it under savepath
        :rtype: Figure object of the created epitope_cluster_plot
    """

    # -----------------------------------------------------------------------------------
    # Create peptide_per_position dictionary which holds all peptides and their position in the protein sequence
    # -----------------------------------------------------------------------------------

    peptide_index = epitopeResult.index.levels[0]
    protein_seq = None
    peptide_per_pos = {}
    for i in peptide_index:
        if transcript_id in i.proteins:
            if protein_seq is None:
                protein_seq = i.proteins[transcript_id]
            for pos in i.get_protein_positions(transcript_id):
                peptide_per_pos[pos] = i

    # -----------------------------------------------------------------------------------
    # Initialize a new line of zeroes
    # -----------------------------------------------------------------------------------

    def zero_add(array, n):
        array.append(np.zeros(n))

    # -----------------------------------------------------------------------------------
    # Create dictionary of alleles with their epitope lists
    # -----------------------------------------------------------------------------------

    ecp = {}
    for allele in epitopeResult.columns:
        ecp_a = []
        previous_pos = -1
        previous_length = -1
        counter = -1
        for pos in xrange(len(protein_seq)):  # iterate over epitopes and their positions
            epitope = peptide_per_pos.get(pos, None)
            if epitope is not None:  # if epitope exists
                epi_value = epitopeResult.loc[(epitope, method), allele]
                if cmp(epi_value, threshold):  # check value vs threshold
                    if previous_pos <= pos < previous_pos+previous_length:  # Am I still within a overlapping region?
                        counter += 1  # jump in next row
                        if len(ecp_a)-1 < counter:  # check if ecp_a needs new row
                            zero_add(ecp_a, len(protein_seq))   # expand ecp_a
                            for epi_pos in xrange(pos, pos+len(epitope)-1):  # write new epitope
                                ecp_a[counter][epi_pos] = max(0, 1 - math.log(epi_value, 50000))
                        else:  # if ecp_a is big enough
                            for epi_pos in xrange(pos, pos+len(epitope)-1):  # write new epitope
                                ecp_a[counter][epi_pos] = max(0, 1 - math.log(epi_value, 50000))
                    previous_pos = pos
                    previous_length = len(epitope)
        if sum(sum(v) for v in ecp_a) > 0:
            ecp[allele] = ecp_a   # insert list of epitopes for new allele in dictionary

    # -----------------------------------------------------------------------------------
    # Define plot
    # -----------------------------------------------------------------------------------

    row_height = 0.25
    nof_rows_per_plot = [len(e) for k, e in sorted(ecp.iteritems())]
    max_rows = max(nof_rows_per_plot)
    nof_rows = sum(nof_rows_per_plot)
    gs = gridspec.GridSpec(len(ecp), 1, height_ratios=[p_height/max_rows for p_height in nof_rows_per_plot])
    fig = plt.figure(figsize=(len(protein_seq)*0.25, row_height*nof_rows))

    cmap = plt.cm.get_cmap('Reds')
    cmap.set_bad(color='0.75', alpha=None)

    column_labels = list(protein_seq)
    x_label_distance = 0.375
    for i, (allele, values) in enumerate(sorted(ecp.iteritems())):
        axs = plt.subplot(gs[i])
        matrix = np.array(values)  # create numpy array with values of dictionary
        if cm_max and cm_min is None:
            cm = axs.pcolormesh(matrix, cmap=cmap, vmax=matrix.max(), vmin=matrix.min())
        else:
            cm = axs.pcolormesh(matrix, cmap=cmap, vmax=cm_max, vmin=cm_min)
        LABEL_SIZE = 10
        axs.xaxis.tick_top()
        axs.set_xticks(np.arange(matrix.shape[1]))
        axs.set_yticks(np.arange(matrix.shape[0]))
        axs.set_yticklabels('')
        axs.set_xticklabels('')
        axs.set_ylim(matrix.shape[0], 0)
        axs.tick_params(axis='x', labelsize=LABEL_SIZE)
        axs.tick_params(axis='y', labelsize=LABEL_SIZE)
        axs.set_xlim(left=0, right=matrix.shape[1])
        axs.text(-3.25, len(matrix)/2, str(allele), size=10)  # write xaxis and yaxis labels as text
        for j in xrange(len(protein_seq)):
            axs.text(j+x_label_distance, -len(matrix)*0.01, str(column_labels[j]), size=10)
        axs.grid(b=True, which='major', axis='both', color='black', linestyle='-')

    cbaxes = fig.add_axes([0.92, 0.683, 0.02, nof_rows_per_plot[0]*0.15])
    cb = plt.colorbar(cm, cax=cbaxes)
    cb.set_ticks([cm_min, cm_max])
    cb.ax.tick_params(labelsize=LABEL_SIZE)

    plt.subplots_adjust(hspace=0.2)
    fig.savefig(savepath, dp8i=75)
    plt.close(fig)

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math


def epitope_density_plot(epitopeResult, transcript_id, method=None, cm_max=None, cm_min=None, show_additional_plot=True,
                         savepath='out.pdf'):
    """

        :param epitopeResult: Pandas Dataframe of an epitope prediction
        as :class:`~Fred2.Core.Result.EpitopePredictionResult` object
        :type epitopeResult: :class:`~Fred2.Core.Base.AEpitopePredictionResult`
        :param str transcript_id: The unique transcript ID of the :class:`~Fred2.Core.Protein.Protein` used in the plot
        :param method: Method used for epitope prediction
        from :class:`~Fred2.EpitopePrediction.PSSM.APSSMEpitopePrediction` resulting
        in an  :class:`~Fred2.Core.Result.EpitopePredictionResult` object
        :type method:
        :param int cm_max: Maximum density for epitope_density_plot, default=None
        :param int cm_min: Minimum density for epitope_density_plot, default=None
        :param boolean show_additional_plot: option to turn on/off additional plot
        :param savepath: Path and name of the file to be saved
        :results: Draws an epitope_density_plot for a given EpitopePredictionResult object and saves it as savepath
        :rtype: Figure object of the created epitope_density_plot
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
            peptide_per_pos[i] = i.get_protein_positions(transcript_id)

    # -----------------------------------------------------------------------------------
    # Initialize and fill epitope_matrix with values
    # -----------------------------------------------------------------------------------

    epitope_matrix = np.zeros((len(epitopeResult.columns), len(protein_seq)))
    for a, allele in enumerate(epitopeResult.columns):
        for key, values in peptide_per_pos.iteritems():
            for i in values:
                for k in xrange(i, i + len(key) - 1):
                    epitope_matrix[a][k] += max(0, 1 - math.log(
                        epitopeResult.loc[(key, epitopeResult.index.levels[1][0]), allele], 50000))

    # -----------------------------------------------------------------------------------
    # Define plot
    # -----------------------------------------------------------------------------------

    cmap = plt.cm.get_cmap('Reds')
    cmap.set_bad(color='0.75', alpha=None)
    allele_num = epitope_matrix.shape[0]
    inch_per_box = 1.25
    fig = plt.figure(figsize=(len(protein_seq), allele_num * inch_per_box))
    gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[allele_num / 3 * inch_per_box, 1],
                                      width_ratios=[len(protein_seq), 1])
    ax1 = plt.subplot(gs[0, :])
    if cm_max and cm_min is None:
        cm = plt.pcolormesh(epitope_matrix, cmap=cmap, vmax=epitope_matrix.max(), vmin=epitope_matrix.min())
    else:
        cm = plt.pcolormesh(epitope_matrix, cmap=cmap, vmax=cm_max, vmin=cm_min)

    LABEL_SIZE = 10
    ax1.xaxis.tick_top()
    ax1.set_xticks(np.arange(epitope_matrix.shape[1]))
    ax1.set_yticks(np.arange(epitope_matrix.shape[0]))
    ax1.set_ylim(bottom=0, top=epitope_matrix.shape[0])
    ax1.set_xlim(left=0, right=epitope_matrix.shape[1])
    ax1.set_yticklabels('')
    ax1.set_xticklabels('')
    plt.tick_params(axis='x', labelsize=LABEL_SIZE)
    plt.tick_params(axis='y', labelsize=LABEL_SIZE)
    y_label_distance = 1
    x_label_distance = 0.275
    column_labels = list(protein_seq)
    for i in range(epitope_matrix.shape[0]):
        # write yaxis labels as text
        ax1.text(-3.5, y_label_distance - inch_per_box / 2 - 0.15, str(epitopeResult.columns.values[i]), size=40)
        y_label_distance += 1
    for i in range(epitope_matrix.shape[1]):
        # write xaxis labels as text
        ax1.text(x_label_distance, allele_num + 0.15, str(column_labels[i]), size=40)
        x_label_distance += 1
    plt.grid(b=True, which='major', axis='both', color='black', linestyle='-')

    top_adjust = 0.9
    cbaxes = fig.add_axes([0.92, top_adjust * 0.250, 0.012, allele_num * 0.0195 + inch_per_box * 0.075])
    if cm_min is None:
        cm_min = epitope_matrix.min()
    if cm_max is None:
        cm_max = epitope_matrix.max()
    cb = plt.colorbar(cm, cax=cbaxes)
    cb.set_ticks([cm_min, cm_max])
    cb.ax.xaxis.set_ticks_position("none")
    cb.ax.yaxis.set_ticks_position("none")
    cb.ax.tick_params(labelsize=40)

    if show_additional_plot is True:
        ax2 = plt.subplot(gs[1, :])
        epitope_matrix_2_y = np.sum(epitope_matrix, axis=0)
        ax2.fill_between(np.arange(len(protein_seq)) + 0.5, epitope_matrix_2_y, facecolor='#F97E60')

    plt.subplots_adjust(top=top_adjust)
    plt.savefig(savepath, dpi=100)
    plt.close(fig)
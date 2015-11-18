from Fred2.Core import Protein, Generator, Allele
import Fred2.Core.Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
import Fred2.Core.Allele
import Fred2.EpitopePrediction.PSSM as PSSM_df
import matplotlib
import numpy as np
import matplotlib.pyplot as plt



def epitope_density_plot(epitopeResult, transcript_id, method=None):
    """
        Description
    """
    # -----------------------------------------------------------------------------------
    # what protein is the peptide from
    # -----------------------------------------------------------------------------------
    s = epitopeResult
    print(s)
    epitope_matrix = np.zeros(shape=(10,10))
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
    ax1.plot([1,2])
    plt.title('Epitope_Density_Plot')
    #plt.show()
    # -----------------------------------------------------------------------------------
    # define epitope matrix
    # -----------------------------------------------------------------------------------
    """
    test_peptide = Fred2.Core.Peptide()
    test_allele = Fred2.Core.Allele()
    matrix_width = len(test_peptide.test_getitem())
    epitope_matrix = PSSM_df.predict(test_peptide, test_allele)
    cm = plt.pcolormesh(epitope_matrix ,cmap=cmap ,vmax=4 ,vmin=0)
    """


protein = Protein("ASDERWQTGHKILPMNVFCY", gene_id=1, transcript_id="someID")
peps = Generator.generate_peptides_from_proteins(protein, 9)
result = EpitopePredictorFactory("BIMAS").predict(peps, alleles=Allele("HLA-A*02:01"))

epitope_density_plot(result, "someID",method="BIMAS")





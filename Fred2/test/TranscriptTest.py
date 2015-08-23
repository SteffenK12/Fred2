from unittest import TestCase

__author__ = 'walzer'

from Fred2.Core import Transcript
from Fred2.Core import Variant
from Fred2.Core import VariationType
from Fred2.Core import MutationSyntax

class TestTranscript(TestCase):
    def setUp(self):
        self.simple = Transcript("")
        self.simple_new = Transcript("")
        self.w_gid = Transcript("", _gene_id="123")
        self.w_tid = Transcript("", _transcript_id="tid")
        self.w_id = Transcript("", "gid", "tid")
        gcg_v1 = Variant("rs5650", VariationType.SNP, 2, 162145588, 'G', 'T',
                         {"NM_002054.4": MutationSyntax("NM_002054.4", 344, 115, "c.344C>A", "p.A115D")}, False, False)

        self.gcg_ts = "gcatagaatgcagatgagcaaagtgagtgggagagggaagtcatttgtaacaaaaactcattatttacagatgagaaatttatattgtcagcgtaatatctgtgaggctaaacagagctggagagtatataaaagcagtgcgccttggtgcagaagtacagagcttaggacacagagcacatcaaaagttcccaaagagggcttgctctctcttcacctgctctgttctacagcacactaccagaagacagcagaaatgaaaagcatttactttgtggctggattatttgtaatgctggtacaaggcagctggcaacgttcccttcaagacacagaggagaaatccagatcattctcagcttcccaggcagacccactcagtgatcctgatcagatgaacgaggacaagcgccattcacagggcacattcaccagtgactacagcaagtatctggactccaggcgtgcccaagattttgtgcagtggttgatgaataccaagaggaacaggaataacattgccaaacgtcacgatgaatttgagagacatgctgaagggacctttaccagtgatgtaagttcttatttggaaggccaagctgccaaggaattcattgcttggctggtgaaaggccgaggaaggcgagatttcccagaagaggtcgccattgttgaagaacttggccgcagacatgctgatggttctttctctgatgagatgaacaccattcttgataatcttgccgccagggactttataaactggttgattcagaccaaaatcactgacaggaaataactatatcactattcaagatcatcttcacaacatcacctgctagccacgtgggatgtttgaaatgttaagtcctgtaaatttaagaggtgtattctgaggccacattgctttgcatgccaataaataaattttcttttagtgttgtgtagccaaaaattacaaatggaataaagttttatcaaaatattgctaaaatatcagctttaaaatatgaaagtgctagattctgttattttcttcttattttggatgaagtaccccaacctgtttacatttagcgataaaattatttttctatgatataatttgtaaatgtaaattattccgatctgacatatctgcattataataataggagaatagaagaactggtagccacagtggtgaaattggaaagagaactttcttcctgaaacctttgtcttaaaaatactcagctttcaatgtatcaaagatacaattaaataaaattttcaagcttctttaccattgtct"
        self.w_v = Transcript(self.gcg_ts, 'GLUC_HUMAN', "NM_002054.4", [gcg_v1])

    def test_consistency(self):
        """
        tests all __*__ (including init)
        test has several asserts! If one fails, the following will not be evaluated!
        """
        self.assertTrue(repr(self.simple) == "TRANSCRIPT: 0\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(self.simple_new) == "TRANSCRIPT: 1\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(self.w_gid) == "TRANSCRIPT: 2\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(self.w_tid) == "TRANSCRIPT: tid\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(self.w_id) == "TRANSCRIPT: tid\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(self.w_v) == "TRANSCRIPT: NM_002054.4\n\tVARIANTS:\n\t\tpos 344: Variant(g.162145588G>T)\n\tSEQUENCE: gcatagaatgcagatgagcaaagtgagtgggagagggaagtcatttgtaacaaaaactcattatttacagatgagaaatttatattgtcagcgtaatatctgtgaggctaaacagagctggagagtatataaaagcagtgcgccttggtgcagaagtacagagcttaggacacagagcacatcaaaagttcccaaagagggcttgctctctcttcacctgctctgttctacagcacactaccagaagacagcagaaatgaaaagcatttactttgtggctggattatttgtaatgctggtacaaggcagctggcaacgttcccttcaagacacagaggagaaatccagatcattctcagcttcccaggcagacccactcagtgatcctgatcagatgaacgaggacaagcgccattcacagggcacattcaccagtgactacagcaagtatctggactccaggcgtgcccaagattttgtgcagtggttgatgaataccaagaggaacaggaataacattgccaaacgtcacgatgaatttgagagacatgctgaagggacctttaccagtgatgtaagttcttatttggaaggccaagctgccaaggaattcattgcttggctggtgaaaggccgaggaaggcgagatttcccagaagaggtcgccattgttgaagaacttggccgcagacatgctgatggttctttctctgatgagatgaacaccattcttgataatcttgccgccagggactttataaactggttgattcagaccaaaatcactgacaggaaataactatatcactattcaagatcatcttcacaacatcacctgctagccacgtgggatgtttgaaatgttaagtcctgtaaatttaagaggtgtattctgaggccacattgctttgcatgccaataaataaattttcttttagtgttgtgtagccaaaaattacaaatggaataaagttttatcaaaatattgctaaaatatcagctttaaaatatgaaagtgctagattctgttattttcttcttattttggatgaagtaccccaacctgtttacatttagcgataaaattatttttctatgatataatttgtaaatgtaaattattccgatctgacatatctgcattataataataggagaatagaagaactggtagccacagtggtgaaattggaaagagaactttcttcctgaaacctttgtcttaaaaatactcagctttcaatgtatcaaagatacaattaaataaaattttcaagcttctttaccattgtct (mRNA)")

    def test_translate(self):
        gcg_var = "A*NADEQSEWEREVICNKNSLFTDEKFILSA*YL*G*TELESI*KQCALVQKYRA*DTEHIKSSQRGLALSSPALFYSTLPEDSRNEKHLLCGWIICNAGTRQLATFPSRHRGEIQIILSFPGRPTQ*S*SDERGQAPFTGHIHQ*LQQVSGLQACPRFCAVVDEYQEEQE*HCQTSR*I*ETC*RDLYQ*CKFLFGRPSCQGIHCLAGERPRKARFPRRGRHC*RTWPQTC*WFFL**DEHHS**SCRQGLYKLVDSDQNH*QEITISLFKIIFTTSPASHVGCLKC*VL*I*EVYSEATLLCMPINKFSFSVV*PKITNGIKFYQNIAKISALKYESARFCYFLLILDEVPQPVYI*R*NYFSMI*FVNVNYSDLTYLHYNNRRIEELVATVVKLERELSS*NLCLKNTQLSMYQRYN*IKFSSFFTIV"
        self.assertTrue(self.w_v.translate() == gcg_var)
        #http://stackoverflow.com/questions/3892218/how-to-test-with-pythons-unittest-that-a-warning-has-been-thrown

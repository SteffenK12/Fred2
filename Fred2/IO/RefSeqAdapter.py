# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: IO.RefSeqAdapter
   :synopsis: DB-Adapter class for RefSeq
.. moduleauthor:: walzer
"""

import logging

from Bio import SeqIO

from Fred2.IO.ADBAdapter import ADBAdapter

#TODO: refactor ... function based on old code
class RefSeqAdapter(ADBAdapter):
    def __init__(self, prot_file=None, prot_vers=None, mrna_file=None, mrna_vers=None):
        self.refseq_prot = self.load(prot_file)
        self.vers_prot = prot_vers
        self.refseq_mrna = self.load(mrna_file)
        self.vers_mrna = mrna_vers

    def load(self, filename):
        refseq_records = dict()
        with open(filename, "rU") as f:
            for record in SeqIO.parse(f, "fasta"):
                ridv = filter(None, record.id.split('|'))[-1]  # NP_001639.1
                rid = ridv.split('.')[0]  # NP_001639
                if rid not in refseq_records:
                    refseq_records[rid] = record
                    refseq_records[rid].dbxrefs.append(ridv)
                    refseq_records[rid].id = rid
                else:
                    print 'claaaash!!'  # TODO no clashes in v.66 but ever?! use logging.warning or something
        return refseq_records

    def get_product_sequence(self, product_refseq):
        """
        fetches product sequence for the given id
        :param product_refseq: given refseq id
        :return: list of dictionaries of the requested sequence, the respective strand and the associated gene name
        """
        if self.refseq_prot:
            if product_refseq in self.refseq_prot:
                return self.refseq_prot[product_refseq]
            else:
                logging.warning('no such sequence')
        else:
            logging.warning('no sequences loaded')

    def get_transcript_sequence(self, transcript_refseq):
        """
        Fetches transcript sequence for the given id
        :param transcript_refseq:
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        if self.refseq_mrna:
            if transcript_refseq in self.refseq_mrna:
                return self.refseq_mrna[transcript_refseq]
            else:
                logging.warning('no such sequence')
        else:
            logging.warning('no sequences loaded')

    def get_transcript_information(self, transcript_refseq):
        pass

# from Bio import Entrez
# >>> rec = Entrez.read(Entrez.esearch(db="protein", term="NP_001005218" ))
# /usr/local/lib/python2.7/dist-packages/biopython-1.64-py2.7-linux-x86_64.egg/Bio/Entrez/__init__.py:451: UserWarning:
# Email address is not specified.
#
# To make use of NCBI's E-utilities, NCBI requires you to specify your
# email address with each request.  As an example, if your email address
# is A.N.Other@example.com, you can specify it as follows:
#    from Bio import Entrez
#    Entrez.email = 'A.N.Other@example.com'
# In case of excessive usage of the E-utilities, NCBI will attempt to contact
# a user at the email address provided before blocking access to the
# E-utilities.
#   E-utilities.""", UserWarning)
# >>> fasta = Entrez.efetch(db="protein", id=rec["IdList"][0], rettype="fasta").read()

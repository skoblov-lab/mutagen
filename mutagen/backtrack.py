"""

Backtrack single amino-acid (AA) substitutions to genomic location in HG38 via
EBI's Protein API. AA at non-zero phased splice-sites are not
supported.

Requirements:
    - biopython
    - requests
    - fn

"""

import json
import operator as op
from typing import NamedTuple, Optional

import requests
from Bio.Data import CodonTable
from fn.iters import group_by

REQUEST = "https://www.ebi.ac.uk/proteins/api/coordinates/location/{prot}:{pos}"

_REVMAP = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

_STDTABLE = CodonTable.unambiguous_dna_by_id[1]
STDCODE = {
    acid:  [codon for codon, _ in pairs] for acid, pairs in
    group_by(op.itemgetter(1), _STDTABLE.forward_table.items()).items()
}
_MITOTABLE = CodonTable.unambiguous_dna_by_id[2]
MITOCODE = {
    acid:  [codon for codon, _ in pairs] for acid, pairs in
    group_by(op.itemgetter(1), _MITOTABLE.forward_table.items()).items()
}


Codon = NamedTuple('Region', [
    ('contig', str), ('transcript', str), ('forward', bool), ('number', int),
    ('start', int), ('stop', int), ('tstart', int), ('tstop', int)
])

_extract_location = op.itemgetter('chromosome', 'ensemblTranslationId',
                                  'geneStart', 'geneEnd')


def locate(uniprot_id: str, position: int) -> Optional[Codon]:
    """
    :param uniprot_id: UniProtKB/SwissProt accession
    :param position: a 0-based amino-acid position
    :return:
    """
    url = REQUEST.format(prot=uniprot_id, pos=position+1)
    response = requests.get(url, headers={"Accept": "application/json"})
    if not response.ok:
        return None
    try:
        locations = json.loads(response.text)['locations'][0]
        contig, transcript, genstart, genstop = _extract_location(locations)
    except (KeyError, IndexError):
        return None
    if abs(genstart-genstop) > 3:
        # the codon is formed by exons flanking a non-zero phased intron
        return None
    is_forward = genstart < genstop
    start, stop = (genstart-1, genstop) if is_forward else (genstop-1, genstart)
    tstart, tstop = position*3, position*3 + 3
    return Codon(
        contig, transcript, is_forward, position, start, stop, tstart, tstop
    )





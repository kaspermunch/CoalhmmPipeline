
from copy import deepcopy
from re import sub

class MafMasker(object):

    def __init__(self, junkchars):
        self.junkchars = junkchars

    def mask(self, maf):
        """
        Masks sequence symbols of MAF inplace
        """
        clone = deepcopy(maf)
        for idx in range(clone.count()):
            clone.alignment[idx][6] = sub(r'[{}]'.format(self.junkchars), 'N', clone.alignment[idx][6])
        return clone
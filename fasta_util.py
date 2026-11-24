"""
ALL OF THE CODE IS HANDWRITTEN BY A211474 :>
"""

from typing import *
from collections import deque

class FASTA:
    ALL_AA = {
        'F','L','S','Y', '*', 'C', 'W', 'P', 'H', 'Q', 'R', 'I', 'M', 'T', 'N', 'K', 'V', 'A', 'D', 'E', 'G' 
    }
    
    __RNA_Rev = {
        'A':'U',
        'T':'A',
        'U':'A',
        'C':'G',
        'G':'C',
    }
    
    __DNA_Rev = {
        'A':'T',
        'T':'A',
        'U':'A',
        'C':'G',
        'G':'C',
    }
    
    __Nuc_AACode_Map = {
        'UUU':'F',
        'UUC':'F',
        'UUA':'L',
        'UUG':'L',
        'UCU':'S',
        'UCC':'S',
        'UCA':'S',
        'UCG':'S',
        'UAU':'Y',
        'UAC':'Y',
        'UAA':'*',
        'UAG':'*',
        'UGU':'C',
        'UGC':'C',
        'UGA':'*',
        'UGG':'W',
        'CUU':'L',
        'CUC':'L',
        'CUA':'L',
        'CUG':'L',
        'CCU':'P',
        'CCC':'P',
        'CCA':'P',
        'CCG':'P',
        'CAU':'H',
        'CAC':'H',
        'CAA':'Q',
        'CAG':'Q',
        'CGU':'R',
        'CGC':'R',
        'CGA':'R',
        'CGG':'R',
        'AUU':'I',
        'AUC':'I',
        'AUA':'I',
        'AUG':'M',
        'ACU':'T',
        'ACC':'T',
        'ACA':'T',
        'ACG':'T',
        'AAU':'N',
        'AAC':'N',
        'AAA':'K',
        'AAG':'K',
        'AGU':'S',
        'AGC':'S',
        'AGA':'R',
        'AGG':'R',
        'GUU':'V',
        'GUC':'V',
        'GUA':'V',
        'GUG':'V',
        'GCU':'A',
        'GCC':'A',
        'GCA':'A',
        'GCG':'A',
        'GAU':'D',
        'GAC':'D',
        'GAA':'E',
        'GAG':'E',
        'GGU':'G',
        'GGC':'G',
        'GGA':'G',
        'GGG':'G',
        
        'TTT':'F',
        'TTC':'F',
        'TTA':'L',
        'TTG':'L',
        'TCT':'S',
        'TCC':'S',
        'TCA':'S',
        'TCG':'S',
        'TAT':'Y',
        'TAC':'Y',
        'TAA':'*',
        'TAG':'*',
        'TGT':'C',
        'TGC':'C',
        'TGA':'*',
        'TGG':'W',
        'CTT':'L',
        'CTC':'L',
        'CTA':'L',
        'CTG':'L',
        'CCT':'P',
        'CAT':'H',
        'CGT':'R',
        'ATT':'I',
        'ATC':'I',
        'ATA':'I',
        'ATG':'M',
        'ACT':'T',
        'AAT':'N',
        'AGT':'S',
        'GTT':'V',
        'GTC':'V',
        'GTA':'V',
        'GTG':'V',
        'GCT':'A',
        'GAT':'D',
        'GGT':'G',
    }
    
    __RNA_Complement_tr = str.maketrans(__RNA_Rev)
    __DNA_Complement_tr = str.maketrans(__DNA_Rev)
    
    
    def __init__(self) -> None:
        self.name:str
        self.seq:str
    
    # @staticmethod
    # def compare2(seq1:str, seq2:str) -> str:
    #     return (''
    #             .join([FASTA.__Nuc_AACode_Map[seq[i-3:i]] for i in range(3 + start,len(seq)+1,3)]))
    
    
    
    @staticmethod
    def toAminoAcid(seq:str, start = 0) -> str:
        return (''
                .join([FASTA.__Nuc_AACode_Map[seq[i-3:i]] for i in range(3 + start,len(seq)+1,3)]))
    
    @staticmethod
    def complementDNA(seq:str, start:int=0) -> str:
        return seq[start::].translate(FASTA.__DNA_Complement_tr)
        
    @staticmethod
    def complementRNA(seq:str, start:int=0) -> str:
        return seq[start::].translate(FASTA.__RNA_Complement_tr)
    
    @staticmethod
    def reverseNA(seq:str) -> str:
        return seq[::-1]
    
    @staticmethod
    def reverseComplementDNA(seq:str, start:int=0) -> str:
        return seq[-start-1:0:-1].translate(FASTA.__DNA_Complement_tr)
    
    @staticmethod
    def reverseComplementRNA(seq:str, start:int=0) -> str:
        return seq[-start-1:0:-1].translate(FASTA.__RNA_Complement_tr)
    
    
    @staticmethod
    def fromFiles(fileDirs:Sequence[str]) -> 'list[FASTA]':
        fastaR = []
        fasta:FASTA|None = None
        # readB = False
        
        for fileDir in fileDirs:
            with open(fileDir, 'r') as text:
                for line in text.readlines():
                    
                    if (line.startswith('>')):
                        fasta = FASTA()
                        fasta.name = line[1:-1]
                        fasta.seq = ""
                        # readB = True
                        fastaR.append(fasta)
                        continue
                        
                    if fasta is not None:
                        fasta.seq += line.strip()
                    
        return fastaR
    
    def read(self):
        return f">{self.name}\n{self.seq}"
    
if __name__ == "__main__":
    import Main
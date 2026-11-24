"""
ALL OF THE CODE IS HANDWRITTEN BY A211474 :>
"""

from fasta_util import *

def Solve(ref_MyGene_FASTA:FASTA, ref_Translated_FASTA:Sequence[FASTA], ref_ANSWERS:list[str]):
    
    """ ============= Question 1 ============= """
    seq_Count = len(ref_MyGene_FASTA.seq)
    ref_ANSWERS[0] = str(seq_Count)
    
    """ ============= Question 2  ============= """
    ATCG_Dict = {}
    
    for codon in ('A','C','G','T'):
        ATCG_Dict[codon] = ref_MyGene_FASTA.seq.count(codon)    
    ref_ANSWERS[1] = '\n'.join([f'{k}: {v}' for k,v in ATCG_Dict.items()])
    
    """ ============= Question 3 ============= """
    ref_ANSWERS[2] = f"{round((ATCG_Dict['G']+ATCG_Dict['C'])/seq_Count*100,2)}%"
    
    """ ============= Question 4 ============= """
    stop_Dict = {fasta.name:fasta.seq.count('*') for fasta in ref_Translated_FASTA}
    ref_ANSWERS[3] = ("Name : StopCodonCount\n")
    ref_ANSWERS[3] += "\n".join(f"  {'>>' if v >= 5 else ''}  {k}:{v}" for k,v in stop_Dict.items())
    ref_ANSWERS[3] += f"\n{len([0 for k,v in stop_Dict.items() if v >= 5])} frames have at least 5 Stop codon"
    
    """ ============= Question 5 ============= """
    seqCount_Dict = {fasta.name:(len(fasta.seq), len(fasta.seq.replace('*',''))) for fasta in ref_Translated_FASTA}
    
    minI = 9999999999
    minN = -1
    
    for i in seqCount_Dict.values():
        diff = abs(i[0] - 446)
        if diff < minI:
            minI = diff
            minN = i[0]
        
    ref_ANSWERS[4] = ("Name : AminoAcidCount inc. * : AminoAcidCount not inc. * \n")
    ref_ANSWERS[4] += "\n".join(f"  {'>>' if v[0] == minN  else ''}  {k}:{v[0]}:{v[1]}" for k,v in seqCount_Dict.items())
        

if __name__ == "__main__":
    import Main
    
    
    


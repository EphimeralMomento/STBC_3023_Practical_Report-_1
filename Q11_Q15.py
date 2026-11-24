"""
ALL OF THE CODE IS HANDWRITTEN BY A211474 :>
"""

from fasta_util import *

def Solve(ref_MyGene_FASTA:FASTA, ref_Translated_FASTA:Sequence[FASTA], ref_ANSWERS:list[str]):
    
    """ ============= Question 11 ============= """
    multipleSegment = ref_Translated_FASTA[0].seq.split('*')
    k = max(multipleSegment, key=lambda x:len(x))
    
    
    ref_ANSWERS[10] = f"> Splitted Foward-Frame +1: {len(k)}"
    for segment in multipleSegment:
        l = len(segment)
        ref_ANSWERS[10] += f"\n  {'>>' if segment==k else ''}  {segment}: {l}"
    
    # ref_ANSWERS[10] += f"\nLongest Segment {k}: {len(k)}"
        
    """ ============= Question 12 ============= """
    ref_ANSWERS[11] = ""
    for revI in range(3):
        multipleSegment = ref_Translated_FASTA[3+revI].seq.split('*')
        k = max(multipleSegment, key=lambda x:len(x))
        
        ref_ANSWERS[11] += f"\n> Splitted Reverse-Complement-Frame -{revI+1}:{len(k)}"
        for segment in multipleSegment:
            l = len(segment)
            ref_ANSWERS[11] += f"\n  {'>>' if segment==k else ''}  {segment}: {l}"
        # ref_ANSWERS[11] += f"\nLongest Segment Reverse-Complement-Frame -{revI+1}: {len(k)}"
        ref_ANSWERS[11] += '\n'
            
            
    ref_ANSWERS[11] = ref_ANSWERS[11][1:-1]
    
    """ ============= Question 13 ============= """
    Found:dict[str, int] = {}
    seq = ref_Translated_FASTA[0].seq
    for c in seq:
        if c in Found:
            Found[c] += 1
        else:
            Found[c] = 1
    
    del(Found['*'])
    m = max(Found.values())
    ref_ANSWERS[12] = "AminoAcid : Count"
    ref_ANSWERS[12] += '\n' + '\n'.join([f'  {'>>' if v == m else ''}  {k}:{v}' for k,v in Found.items()])
    
    """ ============= Question 14 ============= """
    MCount = {fasta:fasta.seq.count('M') for fasta in ref_Translated_FASTA}
    m = max(Found.values())
    ref_ANSWERS[13] = '\n'.join([f'  {'>>' if v > 0 else ''}>  {k.name}:{v}' for k,v in MCount.items()])
    
    """ ============= Question 15 (Related To Question 7?!!) ============= """
    stop_position_1st_pos = 3*34 + 1
    codon3 = ref_MyGene_FASTA.seq[stop_position_1st_pos:stop_position_1st_pos+3]
        
    ref_ANSWERS[14] = f"stop_position_1st_pos: {stop_position_1st_pos+1}"
    ref_ANSWERS[14] += f"\ncodon ({stop_position_1st_pos+1} - {stop_position_1st_pos+3}): ({(codon3[0])}){codon3[1:3]}"
    


    
    
    
if __name__ == "__main__":
    import Main
    
    
    


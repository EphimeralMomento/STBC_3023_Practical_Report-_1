"""
ALL OF THE CODE IS HANDWRITTEN BY A211474 :>
"""

from fasta_util import *

def Solve(ref_MyGene_FASTA:FASTA, ref_Translated_FASTA:Sequence[FASTA], ref_ANSWERS:list[str]):
    
    """ ============= Question 16 ============= """
    ref_MyGene_FASTA_ReSeq = ref_MyGene_FASTA.seq
    ref_MyGene_FASTA_ReSeq = ref_MyGene_FASTA_ReSeq[:-1] + 'C'

    ref_ANSWERS[15] = ""
    for fwdI in range(3):
        actual = ref_Translated_FASTA[fwdI].seq
        regenerated = FASTA.toAminoAcid(ref_MyGene_FASTA_ReSeq, fwdI)
        
        valid = True
        regenDiff = ""
        actualDiff = ""
        for a,r in zip(actual.strip(), regenerated.strip()):
            if (a == r):
                actualDiff += a
                regenDiff += r
            else:
                valid = False
                actualDiff += f"({a})"
                regenDiff += f"({r})"
            
        ref_ANSWERS[15] += f"\n> Regenerated Forward-Frame +{fwdI+1}:"
        ref_ANSWERS[15] += f"\n{regenDiff}"
        ref_ANSWERS[15] += f"\n> translated_frames Forward-Frame +{fwdI+1}:"
        ref_ANSWERS[15] += f"\n{actualDiff}"
        ref_ANSWERS[15] += f"\nValidity: {'Valid' if valid else 'InValid'}\n"
        
    for revI in range(3):
        actual = ref_Translated_FASTA[revI+3].seq
        regenerated = FASTA.toAminoAcid(FASTA.reverseComplementDNA(ref_MyGene_FASTA_ReSeq,revI))
        
        valid = True
        regenDiff = ""
        actualDiff = ""
        for a,r in zip(actual.strip(), regenerated.strip()):
            if (a == r):
                actualDiff += a
                regenDiff += r
            else:
                valid = False
                actualDiff += f"({a})"
                regenDiff += f"({r})"
            
        ref_ANSWERS[15] += f"\n> Regenerated Reverse-Complement-Frame -{revI+1}:"
        ref_ANSWERS[15] += f"\n{regenDiff}"
        ref_ANSWERS[15] += f"\n> translated_frames Reverse-Complement-Frame -{revI+1}:"
        ref_ANSWERS[15] += f"\n{actualDiff}"
        ref_ANSWERS[15] += f"\nValidity: {'Valid' if valid else 'InValid'}\n"
        
    ref_ANSWERS[15] =ref_ANSWERS[15][1:-1]
    
    
    """ ============= Question 17 ============= """
    ref_MyGene_FASTA_ReSeq = ref_MyGene_FASTA.seq
    ref_MyGene_FASTA_ReSeq = 'G' + ref_MyGene_FASTA_ReSeq[1:]
    reTranslated = FASTA.toAminoAcid(ref_MyGene_FASTA_ReSeq)
    ref_ANSWERS[16] = f"First mdoified AminoAcid: {reTranslated[0]}"
    
    
    
    """ ============= Question 18 ============= """
    ref_MyGene_FASTA_ReSeq = ref_MyGene_FASTA.seq
    ref_MyGene_FASTA_ReSeq = ref_MyGene_FASTA_ReSeq[0] + 'C' + ref_MyGene_FASTA_ReSeq[2:]

    ref_ANSWERS[17] = ""
    for fwdI in range(3):
        actual = ref_Translated_FASTA[fwdI].seq
        regenerated = FASTA.toAminoAcid(ref_MyGene_FASTA_ReSeq, fwdI)
        
        valid = True
        regenDiff = ""
        actualDiff = ""
        for a,r in zip(actual.strip(), regenerated.strip()):
            if (a == r):
                actualDiff += a
                regenDiff += r
            else:
                valid = False
                actualDiff += f"({a})"
                regenDiff += f"({r})"
            
        ref_ANSWERS[17] += f"\n> Regenerated Forward-Frame +{fwdI+1}:"
        ref_ANSWERS[17] += f"\n{regenDiff}"
        ref_ANSWERS[17] += f"\n> translated_frames Forward-Frame +{fwdI+1}:"
        ref_ANSWERS[17] += f"\n{actualDiff}"
        ref_ANSWERS[17] += f"\nValidity: {'Valid' if valid else 'InValid'}\n"
        
    for revI in range(3):
        actual = ref_Translated_FASTA[revI+3].seq
        regenerated = FASTA.toAminoAcid(FASTA.reverseComplementDNA(ref_MyGene_FASTA_ReSeq,revI))
        
        valid = True
        regenDiff = ""
        actualDiff = ""
        for a,r in zip(actual.strip(), regenerated.strip()):
            if (a == r):
                actualDiff += a
                regenDiff += r
            else:
                valid = False
                actualDiff += f"({a})"
                regenDiff += f"({r})"
            
        ref_ANSWERS[17] += f"\n> Regenerated Reverse-Complement-Frame -{revI+1}:"
        ref_ANSWERS[17] += f"\n{regenDiff}"
        ref_ANSWERS[17] += f"\n> translated_frames Reverse-Complement-Frame -{revI+1}:"
        ref_ANSWERS[17] += f"\n{actualDiff}"
        ref_ANSWERS[17] += f"\nValidity: {'Valid' if valid else 'InValid'}\n"
        
    ref_ANSWERS[17] =ref_ANSWERS[17][1:-1]
    
    
    """ ============= Question 19 ============= """
    seqCount_Dict = {fasta.name:(len(fasta.seq), len(fasta.seq.replace('*',''))) for fasta in ref_Translated_FASTA}
        
    ref_ANSWERS[18] = ("Name : AminoAcidCount inc. * : AminoAcidCount not inc. * : Diff \n")
    ref_ANSWERS[18] += "\n".join(f"  {'>>' if v[0]-v[1] == 5  else ''}  {k} : {v[0]} : {v[1]} : {v[0]-v[1]}" for k,v in seqCount_Dict.items())
            
    """ ============= Question 20 ============= """
    frm1Dif = ""
    frm2Dif = ""
    
    equal = 0
    diff = 1
    for c1,c2 in zip(ref_Translated_FASTA[0].seq, ref_Translated_FASTA[1].seq):
        if c1 == c2:
            frm1Dif += c1
            frm2Dif += c2
            equal += 1
        else:
            frm1Dif += f'({c1})'
            frm2Dif += f'({c2})'
            diff += 1
        
    ref_ANSWERS[19] = f"Frame +1: {frm1Dif}"
    ref_ANSWERS[19] += f"\nFrame +2: {frm2Dif}"
    ref_ANSWERS[19] += f"\nEquals: {equal}"
    ref_ANSWERS[19] += f"\nDiffs.: {diff}"
                
        
        
    
    
    
if __name__ == "__main__":
    import Main
    
    
    



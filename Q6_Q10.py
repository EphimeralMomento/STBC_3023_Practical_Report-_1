from fasta_util import *

def Solve(ref_MyGene_FASTA:FASTA, ref_Translated_FASTA:Sequence[FASTA], ref_ANSWERS:list[str]):
    
    """ ============= Question 6 ============= """
    rev3 = FASTA.toAminoAcid(FASTA.reverseComplementRNA(ref_MyGene_FASTA.seq, 2))
    ref_ANSWERS[5] = f"> Regenerated Reverse-Complement-Frame -3:"
    ref_ANSWERS[5] += f'\n{rev3}'
    ref_ANSWERS[5] += f'\nTryptophan Count: {rev3.count('W')}'

    """ ============= Question 7 ============= """
    fwd2 = FASTA.toAminoAcid(ref_MyGene_FASTA.seq[1::])
    fwd2AA = fwd2.find('*')
    stop_position_1st_pos = 3*fwd2.find('*') + 1
    codon3 = ref_MyGene_FASTA.seq[stop_position_1st_pos:stop_position_1st_pos+3]
    
    ref_ANSWERS[6] = f"stop_position_aminoAcid: {fwd2AA+1}"
    ref_ANSWERS[6] += f"\nstop_position_1st_pos: {stop_position_1st_pos+1}"
    ref_ANSWERS[6] += f"\ncodon: ({(codon3[0])}){codon3[1:3]}"
    
    """ ============= Question 8 ============= """
    stop_Dict = {fasta:fasta.seq.count('*') for fasta in ref_Translated_FASTA}
    ref_ANSWERS[7] = ("Name : StopCodonCount\n")
    ref_ANSWERS[7] = "\n".join(f"{'>>' if v >= 15 else ''} {k.name}:{v}" for k,v in stop_Dict.items())
    
    """ ============= Question 9 (ALL VALID??) ============= """
    ref_ANSWERS[8] = ""
    for fwdI in range(3):
        actual = ref_Translated_FASTA[fwdI].seq
        regenerated = FASTA.toAminoAcid(ref_MyGene_FASTA.seq, fwdI)
        
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
            
        ref_ANSWERS[8] += f"\n> Regenerated Forward-Frame +{fwdI+1}:"
        ref_ANSWERS[8] += f"\n{regenDiff}"
        ref_ANSWERS[8] += f"\n> translated_frames Forward-Frame +{fwdI+1}:"
        ref_ANSWERS[8] += f"\n{actualDiff}"
        ref_ANSWERS[8] += f"\nValidity: {'Valid' if valid else 'InValid'}\n"
        
    for revI in range(3):
        actual = ref_Translated_FASTA[revI+3].seq
        regenerated = FASTA.toAminoAcid(FASTA.reverseComplementDNA(ref_MyGene_FASTA.seq,revI))
        
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
            
        ref_ANSWERS[8] += f"\n> Regenerated Reverse-Complement-Frame -{revI+1}:"
        ref_ANSWERS[8] += f"\n{regenDiff}"
        ref_ANSWERS[8] += f"\n> translated_frames Reverse-Complement-Frame -{revI+1}:"
        ref_ANSWERS[8] += f"\n{actualDiff}"
        ref_ANSWERS[8] += f"\nValidity: {'Valid' if valid else 'InValid'}\n"
        
    ref_ANSWERS[8] = ref_ANSWERS[8][1:-1]
    
    
    """ ============= Question 10 ============= """
    minI = 9999999999
    minN:FASTA
    
    for k,v in stop_Dict.items():
        if v < minI:
            minI = v
            minN = k
        
    ref_ANSWERS[9] = ("Name : StopCodonCount \n")
    ref_ANSWERS[9] = "\n".join(f"  {'>>' if k == minN  else ''}>  {k.name}:{v}" for k,v in stop_Dict.items())
        

    


    
    
    
if __name__ == "__main__":
    import Main
    
    
    


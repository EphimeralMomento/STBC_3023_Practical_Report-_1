"""
ALL OF THE CODE IS HANDWRITTEN BY A211474 :>
"""



from typing import *
from fasta_util import *
import Q1_Q5
import Q6_Q10
import Q11_Q15
import Q16_Q20


QUESTIONS = [
    "Q01. After parsing My_gene.fasta and removing its FASTA header and line breaks, how many DNA nucleotides are present?",
    "Q02. Write code to count A, C, G and T in My_gene.fasta. Select ALL statements that match the program output.",
    "Q03. What GC percentage does a program calculate for the input DNA, rounded to two decimal places?",
    "Q04. After parsing translated_frames.fasta, how many sequence records have AT LEAST FIVE * stop symbols?",
    "Q05. Which of these reading frames contain exactly 446 amino-acid symbols in translated_frames.fasta, INCLUDING any * stop symbol? Select ALL.",
    "Q06. Regenerate reverse-complement frame -3 using the standard genetic code, then count the number of tryptophan (W) symbols. What is the count?",
    "Q07. What is the 1-based amino-acid position of the FIRST * stop symbol in translated frame +2?",
    "Q08. Which frames contain AT LEAST FIFTEEN * stop symbols in translated_frames.fasta? Select ALL.",
    "Q09. Regenerate all six translations from My_gene.fasta and perform exact string comparisons to the uploaded outputs. Which of these comparison tests return TRUE? Select ALL.",
    "Q10. Which translated reading frame has the FEWEST * stop symbols when each full record is counted?",
    "Q11. For translated frame +1, what is the length of the LONGEST continuous segment without *? Exclude * from the segment length.",
    "Q12. Among the THREE reverse-strand translations (-1, -2, -3), which has the longest continuous stop-free segment?",
    "Q13. Which amino acid is MOST frequent in translated frame +1 (ignore the terminal *), and how many times does it appear?",
    "Q14. Which output frames contain AT LEAST ONE methionine symbol M, anywhere in the translated sequence? Select ALL.",
    "Q15. The first stop in frame +2 occurs at translated residue 35. Which 1-based DNA coordinate interval encodes that STOP triplet in the ORIGINAL forward DNA?",
    "Q16. Make ONE copy of the DNA and mutate the LAST nucleotide (position 1338) from A to C. Regenerate the six frames. Which observations are TRUE? Select ALL.",
    "Q17. Independently reset to the ORIGINAL DNA. Replace nucleotide 1 from A to G (ATG becomes GTG). What is the NEW first amino-acid symbol in frame +1?",
    "Q18. Independently reset to the ORIGINAL DNA. Replace nucleotide 2 from T to C. Compare each complete translated string before and after the change. Which frames DIFFER? Select ALL.",
    "Q19. Create a new protein-file version by REMOVING all * characters from each translated record. Which frames are then EXACTLY FIVE symbols shorter? Select ALL.",
    "Q20. Without introducing gaps or doing a sequence alignment, compare frame +1 and frame +2 character-by-character over their first 445 amino-acid positions. At how many positions are the symbols identical?"
    ]
ANSWERS:list[str] = ["No Answers!"]*20


MY_GENE_DIR = "My_gene.fasta"
TRANSLATED_FRAMES_DIR = "translated_frames.fasta"



My_Gene_fasta = FASTA.fromFiles([MY_GENE_DIR])
Translated_Frames_fasta = FASTA.fromFiles([TRANSLATED_FRAMES_DIR])

Q1_Q5.Solve(My_Gene_fasta[0], Translated_Frames_fasta, ANSWERS)
Q6_Q10.Solve(My_Gene_fasta[0], Translated_Frames_fasta, ANSWERS)
Q11_Q15.Solve(My_Gene_fasta[0], Translated_Frames_fasta, ANSWERS)
Q16_Q20.Solve(My_Gene_fasta[0], Translated_Frames_fasta, ANSWERS)


with open("run_log.txt", 'w') as writeOut:
    for i in range(len(QUESTIONS)):
        res = "\t" + "\n\t".join(ANSWERS[i].split('\n'))
        print(QUESTIONS[i])
        print(res)
        print()
        
        writeOut.write(QUESTIONS[i] + '\n')
        writeOut.write(res + '\n')
        writeOut.write('\n')
        

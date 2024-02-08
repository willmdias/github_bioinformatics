pacman::p_load(pacman, msa, tidyr, dplyr, genepop, PopGenome, BiocManager, 
               Biostrings, seqinr, phangorn)

#7. ####
seq1 <- readDNAStringSet("cytB_Aparallactus-capensis.fasta")
seq1

seq2 <- readDNAStringSet("cytB_Bitis-nasicornis.fasta")
seq2

seq3 <- readDNAStringSet("cytB_Dispholidus-typus.fasta")
seq3

seq4 <- readDNAStringSet("cytB_Lytorhynchus-diadema.fasta")
seq4

seq5 <- readDNAStringSet("cytB_Naja-annulata.fasta")
seq5

#8. ####
sequences <- c(seq1, seq2, seq3, seq4, seq5)
sequences

names(sequences) <- c ("A. capensis", "B. nasicornis", "D. typus", "L. diadema",
                       "N. annulata")
sequences

#9. ####
msaClustalW(sequences)

msaMuscle(sequences)

#10. ####
FirstAlignment <- msa(sequences)
FirstAlignment

print(FirstAlignment, show="complete")

#11. ####
ncol(FirstAlignment)

#12. ####
GC_content <- alphabetFrequency(FirstAlignment)
GC_content

#13. ####
SeqinrAlignment <- msaConvert(FirstAlignment, type="seqinr::alignment")
SeqinrAlignment

#14. ####
DMatrix <- dist.alignment(SeqinrAlignment, "identity")
DMatrix

#15. ####
Protein1 <- Biostrings::translate(seq1, genetic.code=GENETIC_CODE, 
                      no.init.codon=FALSE, if.fuzzy.codon="error")
names(Protein1) <- "A. capensis"
Protein1

#16. ####
Alignment_phyDat <- msaConvert(FirstAlignment, type="phangorn::phyDat")
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")
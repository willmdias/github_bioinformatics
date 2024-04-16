pacman::p_load(Biostrings, msa, evobiR)

ITS <- readDNAStringSet("R_Input/ITS.fasta")
ITS_A <- msaMuscle(ITS)
writeXStringSet(unmasked(ITS_A), "R_Output/ITS_Aligned.fasta", append = FALSE)

rpoB <- readDNAStringSet("R_Input/rpoB.fasta")
rpoB_A <- msaMuscle(rpoB)
writeXStringSet(unmasked(rpoB_A), "R_Output/rpoB_Aligned.fasta", append = FALSE)

trnC <- readDNAStringSet("R_Input/trnC.fasta")
trnC_A <- msaMuscle(trnC)
writeXStringSet(unmasked(trnC_A), "R_Output/trnC_Aligned.fasta", append = FALSE)

trnE <- readDNAStringSet("R_Input/trnE.fasta")
trnE_A <- msaMuscle(trnE)
writeXStringSet(unmasked(trnE_A), "R_Output/trnE_Aligned.fasta", append = FALSE)

trnT <- readDNAStringSet("R_Input/trnT.fasta")
trnT_A <- msaMuscle(trnT)
writeXStringSet(unmasked(trnT_A), "R_Output/trnT_Aligned.fasta", append = FALSE)

setwd("R_Output")
SuperMatrix(missing = "-", prefix = "Concatenated", save = T)

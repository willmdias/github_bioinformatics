# Load Packages ####
pacman::p_load(pacman, msa, tidyr, dplyr, genepop, PopGenome, BiocManager, 
               Biostrings, seqinr, phangorn, UniprotR, protti, 
               GenomicAlignments, r3dmol)

# Set Directory ####
#It seems like using an R project file automatically sets the working directory 
#to whatever folder the project file is located in, removing the need to
#manually set the directory using setwd(). Please correct me if I am wrong. 
# Yup, you're correct

#1. Import and align your DNA sequences. ####
Seqs <- readDNAStringSet("Input/sequences.fasta")
#Reads the .fasta file containing the DNA sequences and saves it as a variable.
Seqs

names(Seqs) <- c ("Individual 1", "Individual 2", "Individual 3", "Individual 4", 
                  "Individual 5", "Individual 6", "Individual 7", "Individual 8", 
                  "Individual 9", "Individual 10", "Individual 11", "Individual 12", 
                  "Individual 13", "Individual 14", "Individual 15", "Individual 16", 
                  "Individual 17", "Individual 18", "Individual 19", "Individual 20")
#Renames the sequences.
# renaming not fully necessary, as the old names are just as long
# there is also the risk of missing up the identifier numbers when adding the new ones manually
Seqs

Alignment <- msaClustalW(Seqs)
#Creates a multiple sequence alignment from the sequence data utilizing the
#ClustalW method.
print(Alignment, show="complete")
#Prints the entire alignment.

#2. Check to see how different your samples are from one another. Are any of 
#them different from the rest? If so, what kinds of mutations do you observe in 
#this individual (or individuals)? ####

#Individual #4 (Homo_sapiens_4) possesses a point mutation (C to A) at bp 39.
#Individual #6 (Homo_sapiens_6) possesses a deletion (of the C that is normally 
#at bp 3) and multiple point mutations (A to G at bp 47, T to C at bp 134,
#A to T at bp 145, C to G at bp 152, G to C at bp 586, and A to G at bp 623).
#Individual #10 (Homo_sapiens_10) possesses two point mutations (C to G at bp 39 
#and A to T at bp 45).

#3.	You suspect that an individual (or individuals) in this population might 
#have some mutations in this gene, but you don’t know what this gene might be. 
#Compare your sequences to a database to figure out what the gene is. Export 
#your data, paste it into the relevant database search engine, and add your 
#results to a comment line in R. What is the gene? What is the accession number 
#of the best match to your search? ####

#The gene is the hemoglobin subunit beta (HBB) gene. The best match for all
#individuals except #6 was "Homo sapiens hbb gene for beta globin, partial cds, 
#note: HbLimassol Cd8(AAG>AAC)" (Accession number LC121775.1). The best match
#for individual #6 was "Homo sapiens mutant hemoglobin beta chain (HBB) gene, 
#partial cds" (Accession number AY356351.1). The search engine used was NCBI's
#BLAST.

#4.	Find the individual that is the most different from the rest of the 
#individuals in your dataset. Translate that sequence to protein. Write it to a 
#fasta file. ####

#As previously stated in questions #2 and #3, Individual #6 was the most 
#different by a considerable margin.
Protein <- Biostrings::translate(Seqs$"Individual 6", genetic.code=GENETIC_CODE, 
                                  no.init.codon=TRUE, if.fuzzy.codon="error")
#Translates the DNA sequence of Individual #6 into a protein and saves it as a
#separate variable. "no.init.codon=TRUE" as the sequence does NOT start with 
#ATG, the start codon.
Protein

ProteinSet <- AAStringSet(Protein)
#Transforms the AAString protein sequence into an AAStringSet. This is needed
#for the writeXStringSet function to work properly.
ProteinSet

writeXStringSet(ProteinSet, "Output/Protein_mutated.fasta", append = FALSE)
#Creates a .fasta file containing the protein sequence and saves it to the 
#"Output" folder.

#5.	Use a database to figure out what your protein matches to. Click on the 
#record for the best match. What is the accession number of this entry? ####

#The best match was a human hemoglobin subunit beta protein (Accession number 
#A0A0J9YWK4). The search engine/database used was UniProt BLAST. This accession
#number was manually saved to a .csv file in the "Input" folder.

#6.	Either using R or by searching in the database, what disease(s) is this gene
#associated with? Does this person have the disease? ####
Accession <- readLines("Input/UniProt1Mutated.csv")
#Reads the accession number in the .csv file and saves it as a character string.
Accession

GO <- GetProteinGOInfo(Accession)
#Fetches gene ontology terms for the protein based on its accession number.
GO
#There are no gene ontology terms associated with this specific protein.
Path <- GetPathology_Biotech(Accession)
#Returns information on any pathologies associated with the protein.
Path
#There are no pathologies associated with this specific protein.
Dis <- Get.diseases(Path)
#Returns information on any diseases associated with the protein.
Dis
#There are no diseases associated with this specific protein.

#Trying again with an example of the protein in its "normal" state (I used 
#Homo_sapiens_1/Individual #1, although any individuals other than #4, #6, 
#and #10 could be used for this since they share the exact same sequence).
Protein2 <- Biostrings::translate(Seqs$"Individual 1", genetic.code=GENETIC_CODE, 
                                 no.init.codon=TRUE, if.fuzzy.codon="error")
#Translates the DNA sequence of Individual #1 into a protein and saves it as a
#separate variable. "no.init.codon=TRUE" as the sequence does NOT start with 
#ATG, the start codon.
ProteinSet2 <- AAStringSet(Protein2)
#Transforms the AAString protein sequence into an AAStringSet. This is needed
#for the writeXStringSet function to work properly.
ProteinSet2

writeXStringSet(ProteinSet2, "Output/Protein_normal.fasta", append = FALSE)
#Creates a .fasta file containing the protein sequence and saves it to the 
#"Output" folder.
#Again using UniProt BLAST, the greatest match (for a human HBB) was A0A2R8Y7R2.

Accession2 <- readLines("Input/UniProt2Normal.csv")
#Reads the accession number in the .csv file and saves it as a character string.
Accession2

GO2 <- GetProteinGOInfo(Accession2)
#Fetches gene ontology terms for the protein based on its accession number.
GO2
#There are no gene ontology terms associated with this specific protein.
Path2 <- GetPathology_Biotech(Accession2)
#Returns information on any pathologies associated with the protein.
Path2
#There are no pathologies associated with this specific protein.
Dis2 <- Get.diseases(Path2)
#Returns information on any diseases associated with the protein.
Dis2
#There are no diseases associated with this specific protein.

#Both the mutated and non-mutated versions of the protein provided in the 
#sequences file did not yield any results. However, it is known that the protein
#in question is a human hemoglobin subunit beta. Therefore, I will be using the 
#best UniProt match for the simple search "HBB Human" (accession no. P68871).
Accession3 <- readLines("Input/UniProt3Example.csv")
#Reads the accession number in the .csv file and saves it as a character string.
Accession3

GO3 <- GetProteinGOInfo(Accession3)
PlotGOAll(GOObj = GO3, Top = 10, directorypath = "Output", width = 8, height = 5)
#Creates a graph from the results of the gene ontology term search.
#Automatically saved to the "Output" folder.

Path3 <- GetPathology_Biotech(Accession3)
#Returns information on any pathologies associated with the protein.
#I did not print the information this time because Get.diseases summarizes the 
#same information in a much more concise manner without creating a massive wall
#of text.
Dis3 <- Get.diseases(Path3)
#Returns information on any diseases associated with the protein.
Dis3
#This protein is associated with Heinz body anemias, Beta-thalassemia, and 
#sickle cell disease.
#Since Individual 6's protein is heavily mutated compared to the rest of the
#individuals in the "sequences.fasta" file, it can be inferred that they have
#one of these diseases.

#7.	What is the 3-dimensional structure of this protein? You can include a 
#screenshot or download of a photo of this structure in your GitHub repository. ####

#MUTATED PROTEIN:
Fetch <- fetch_uniprot(Accession)
#Fetches information about the protein through UniProt.
print(Fetch, show = "complete")
#The protein is not present in the Protein Database (xref_pdb = NA). r3dmol 
#requires either a Protein Database ID or a .pdb file in order to work properly.
#As such, I utilized its .pdb files from AlphaFold. All structures were manually
#saved to the "Output" folder.

Predict <- fetch_alphafold_prediction(Accession)
#Fetches information regarding the 3D structures of the protein.
Predict

#r3dmol creates visualizations of the 3D structures directly in RStudio.
r3dmol(lowerZoomLimit = 1, upperZoomLimit = 350) %>%
  #Sets up the 3D viewer.
  m_add_model(data = "Input/AF-A0A0J9YWK4-F1-model_v4.pdb") %>%
  #Adds the protein to the viewer.
  m_zoom_to() %>%
  #Sets up the zoom level so that the entire protein is visible.
  m_set_style(style = m_style_cartoon(color = "#b5e61d")) %>%
  #Sets a style for the protein.
  m_set_style(sel = m_sel(ss = "h"), style = m_style_cartoon(color = "#475ecc")) %>%
  #Sets a separate style for the alpha helices.
  m_set_style(sel = m_sel(ss = "s"), style = m_style_cartoon(color = "#ed1c24", arrows = FALSE)) %>%
  #Sets a separate style for the beta sheets. The "arrows" function controls the 
  #presence or absence of arrows that indicate in which directions the sheets run.
  m_rotate(angle = 60, axis = "x")
#Rotates the protein by a given angle on a given axis.

#NORMAL PROTEIN:
Fetch2 <- fetch_uniprot(Accession2)
#Fetches information about the protein through UniProt.
print(Fetch2, show = "complete")
#This protein is also not present in the Protein Database (xref_pdb = NA). The 
#.pdb file from AlphaFold will be used again.

Predict2 <- fetch_alphafold_prediction(Accession2)
#Fetches information regarding the 3D structures of the protein.
Predict2

#r3dmol creates visualizations of the 3D structures directly in RStudio.
r3dmol(lowerZoomLimit = 1, upperZoomLimit = 350) %>%
  #Sets up the 3D viewer.
  m_add_model(data = "Input/AF-A0A2R8Y7R2-F1-model_v4.pdb") %>%
  #Adds the protein to the viewer.
  m_zoom_to() %>%
  #Sets up the zoom level so that the entire protein is visible.
  m_set_style(style = m_style_cartoon(color = "#b5e61d")) %>%
  #Sets a style for the protein.
  m_set_style(sel = m_sel(ss = "h"), style = m_style_cartoon(color = "#475ecc")) %>%
  #Sets a separate style for the alpha helices.
  m_set_style(sel = m_sel(ss = "s"), style = m_style_cartoon(color = "#ed1c24", arrows = FALSE)) %>%
  #Sets a separate style for the beta sheets. The "arrows" function controls the 
  #presence or absence of arrows that indicate in which directions the sheets run.
  m_rotate(angle = 60, axis = "x")
#Rotates the protein by a given angle on a given axis.

#P.S.: Please let me know if I have to explain repeated functions that were 
#previously used and explained. I did it just in case, but I felt it was a bit
#redundant to do so.

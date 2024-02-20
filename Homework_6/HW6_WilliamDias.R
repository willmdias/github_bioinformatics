# Load Packages ####
pacman::p_load(pacman, msa, tidyr, dplyr, genepop, PopGenome, BiocManager, 
               Biostrings, seqinr, phangorn, UniprotR, protti, 
               GenomicAlignments, r3dmol)
#This is how I was taught to do this in my Biological Statistics class. If you'd 
#prefer me to use library() for all packages separately, please let me know.
# pacman works great

# Set Directory ####
#It seems like using an R project file automatically sets the working directory 
#to whatever folder the project file is located in, removing the need to
#manually set the directory using setwd(). Please correct me if I am wrong. 

#3. ####
seq <- readDNAStringSet("cytB_Naja-annulata2.fasta")
#Reads the DNA sequence and saves it as a variable.

Protein <- Biostrings::translate(seq, genetic.code = GENETIC_CODE,
                                 no.init.codon = FALSE, if.fuzzy.codon = "error")
#Translates the DNA sequence into a protein and saves it as a separate variable.
#"no.init.codon" was used as the DNA sequence starts with ATG, the start codon.

names(Protein) <- "Naja annulata"
#Renames the sequence, replacing the default GenBank name.
Protein

writeXStringSet(Protein, "Output/Naja-annulata_protein.fasta", append = FALSE)
#Creates a .fasta file containing the protein sequence and saves it to the 
#"Output" folder.

#4. ####
accession <- read.csv("accession.csv", header = FALSE)
# First value got set as the header name, so if we read it in with header = FALSE
# it automatically assigns a header name and correctly reads in all the values
#Reads the .csv file containing the accession numbers of the 5 most similar 
#sequences obtained through a UniProt BLAST search and saves it as a variable.
accession$V1 # adding the $V1 does the same thing as readLines
#P.S.: The 5th most similar protein found through UniProt BLAST (Q9MLJ3) was 
#replaced by the 6th most similar (Q9MLK0) as the fetch_alphafold_prediction 
#function in question #13 returns no information regarding Q9MLJ3. There is no 
#predicted structure for Q9MLJ3 on the AlphaFold website either.

#5. ####
dat <- readLines("accession.csv")
#Saves the data in the .csv file as a character string.
dat

#6. ####
GO <- GetProteinGOInfo(dat)
#Fetches gene ontology terms for the 5 proteins originally from the .csv file.
GO

#7. ####
PlotGoInfo(GO, directorypath = "PlotGoInfo")
#Creates a graph from the results of the gene ontology term search. 
#P.S.: The "Plots" window needs to be significantly enlarged in order for the 
#plot to look correct.
#Manually saved to the "Output" folder.
# adding the directorypath automatically saves the file

#8. ####
PlotGOAll(GOObj = GO, Top = 10, directorypath = "Output", width = 8, height = 5)
#Different way of plotting the gene ontology information - "Handy visualization
#for publications". Automatically saved to the "Output" folder.

#9. What are some interesting GO terms for your gene? ####
#Respiratory electron transport chain, metal ion binding, ubiquinol-cytochrome-c
#reductase activity, mitochondrial inner membrane, and respiratory chain complex
#III.

#10. ####
path <- GetPathology_Biotech(dat)
#Returns information on any pathologies associated with the 5 proteins.
path
#There are no pathologies associated with these proteins.

dis <- Get.diseases(path)
#Returns information on any diseases associated with the 5 proteins.
dis
#No diseases are associated with these proteins. This was also observed through
#the "GetPathology_Biotech()" function ("Involvement.in.disease = NA" for all).

#11. ####
fetch <- fetch_uniprot(dat)
#Fetches information about the 5 proteins through UniProt.
print(fetch, show = "complete")
#These proteins are not present in the Protein Database (xref_pdb = NA). As
#such, the provided "backup" IDs will be used for #12.

#12. ####
backup <- readLines("q12_backupaccession.csv")
#Reads the .csv file containing the 2 "backup" sequences and saves it as a 
#variable.
fetch_pdb(backup)
#Fetches information about the 2 proteins through the Protein Database.

#13. ####
predict <- fetch_alphafold_prediction(dat)
#Fetches information regarding the 3D structures of the proteins.
predict

#r3dmol creates visualizations of the 3D structures directly in RStudio.
#Q9MLJ9:
r3dmol(lowerZoomLimit = 1, upperZoomLimit = 350) %>%
#Sets up the 3D viewer.
  m_add_model(data = "AF-Q9MLJ9-F1-model_v4.pdb") %>%
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
  m_rotate(angle = 90, axis = "y")
#Rotates the protein by a given angle on a given axis.

#r3dmol requires either a Protein Database ID or a .pdb file in order to work
#properly. Since my proteins did not have a Protein Database ID (only present in
#the UniProt database), I utilized their .pdb files from AlphaFold.
#Structures were manually saved to the "Output" folder.

#Q9MLJ1:
r3dmol(lowerZoomLimit = 1, upperZoomLimit = 350) %>%
  m_add_model(data = "AF-Q9MLJ1-F1-model_v4.pdb") %>%
  m_zoom_to() %>%
  m_set_style(style = m_style_cartoon(color = "#b5e61d")) %>%
  m_set_style(sel = m_sel(ss = "h"), style = m_style_cartoon(color = "#475ecc")) %>%
  m_set_style(sel = m_sel(ss = "s"), style = m_style_cartoon(color = "#ed1c24", arrows = FALSE)) %>%
  m_rotate(angle = 90, axis = "y")

#Q9MLK1:
r3dmol(lowerZoomLimit = 1, upperZoomLimit = 350) %>%
  m_add_model(data = "AF-Q9MLK1-F1-model_v4.pdb") %>%
  m_zoom_to() %>%
  m_set_style(style = m_style_cartoon(color = "#b5e61d")) %>%
  m_set_style(sel = m_sel(ss = "h"), style = m_style_cartoon(color = "#475ecc")) %>%
  m_set_style(sel = m_sel(ss = "s"), style = m_style_cartoon(color = "#ed1c24", arrows = FALSE)) %>%
  m_rotate(angle = 90, axis = "y")

#Q9MLK7:
r3dmol(lowerZoomLimit = 1, upperZoomLimit = 350) %>%
  m_add_model(data = "AF-Q9MLK7-F1-model_v4.pdb") %>%
  m_zoom_to() %>%
  m_set_style(style = m_style_cartoon(color = "#b5e61d")) %>%
  m_set_style(sel = m_sel(ss = "h"), style = m_style_cartoon(color = "#475ecc")) %>%
  m_set_style(sel = m_sel(ss = "s"), style = m_style_cartoon(color = "#ed1c24", arrows = FALSE)) %>%
  m_rotate(angle = 90, axis = "y")

#Q9MLK0:
r3dmol(lowerZoomLimit = 1, upperZoomLimit = 350) %>%
  m_add_model(data = "AF-Q9MLK0-F1-model_v4.pdb") %>%
  m_zoom_to() %>%
  m_set_style(style = m_style_cartoon(color = "#b5e61d")) %>%
  m_set_style(sel = m_sel(ss = "h"), style = m_style_cartoon(color = "#475ecc")) %>%
  m_set_style(sel = m_sel(ss = "s"), style = m_style_cartoon(color = "#ed1c24", arrows = FALSE)) %>%
  m_rotate(angle = 90, axis = "y")

#For proteins that do have a Protein Database ID (such as the two provided 
#"backups"), r3dmol works without the need for downloading any external files.

#1ZMR:
r3dmol(lowerZoomLimit = 1, upperZoomLimit = 350) %>%
  m_add_model(data = m_fetch_pdb("1zmr")) %>%
  m_zoom_to() %>%
  m_set_style(style = m_style_cartoon(color = "#b5e61d")) %>%
  m_set_style(sel = m_sel(ss = "h"), style = m_style_cartoon(color = "#475ecc")) %>%
  m_set_style(sel = m_sel(ss = "s"), style = m_style_cartoon(color = "#ed1c24", arrows = FALSE)) %>%
  m_rotate(angle = 90, axis = "y")

#2HWG:
r3dmol(lowerZoomLimit = 1, upperZoomLimit = 350) %>%
  m_add_model(data = m_fetch_pdb("2hwg")) %>%
  m_zoom_to() %>%
  m_set_style(style = m_style_cartoon(color = "#b5e61d")) %>%
  m_set_style(sel = m_sel(ss = "h"), style = m_style_cartoon(color = "#475ecc")) %>%
  m_set_style(sel = m_sel(ss = "s"), style = m_style_cartoon(color = "#ed1c24", arrows = FALSE)) %>%
  m_rotate(angle = 90, axis = "y")
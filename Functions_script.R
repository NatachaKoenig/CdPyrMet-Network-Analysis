#####################################################################
## PROJECT : ProteoGam                                             ##
## STUDIES : Tox network analysis of raw proteomic data            ##
## AUTHOR : Natacha Koenig                                         ##
## DATE : July 2020                                                ##
## SCRIPT : Functions script                                       ##
#####################################################################

#-------------------------------------------------------------------
#  INTRODUCTORY NOTE                            
#-------------------------------------------------------------------
# This script contains the defintion of the fonctions of the analysis
# This script is executed with the command > source("script.R")

#-------------------------------------------------------------------
#  INDEX.SUMMARY OF THE SCRIPT                   
#-------------------------------------------------------------------

# 1.  FUNCTION ADD_PROT_INFOS    

#-------------------------------------------------------------------
#  1. FUNCTION ADD_PROT_INFOS                
#-------------------------------------------------------------------
## ADD COLUMN MODULE COLOR, kVAL AND ME VAL TO ANNOTATION 

#annot_module : the module to annotate
#color : the color of the module
#kval : kTotal, kWithin, kOut, kDiff
#MEval : the ModuleMembership values

ADD_PROT_INFOS <- function(annot_module, color, kval, MEval){
  d=dim(annot_module)
  Module_color=rep(color, d[1])
  annot=data.frame(annot_module, Module=Module_color)
  annot=left_join(annot, kval, by="Contig")
  annot=left_join(annot, MEval, by="Contig")
  return(annot)
}


print("All the functions have been executed !")
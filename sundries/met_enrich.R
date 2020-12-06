mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-c("Acetoacetic acid","Beta-Alanine","Creatine","Dimethylglycine","Fumaric acid","Glycine","Homocysteine","L-Cysteine","L-Isolucine","L-Phenylalanine","L-Serine","L-Threonine","L-Tyrosine","L-Valine","Phenylpyruvic acid","Propionic acid","Pyruvic acid","Sarcosine")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PerformDetailMatch(mSet, "L-Isolucine");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "L-Isolucine", "L-Isoleucine");

## S/N filter
a<-read.csv("G:/test/LipidSearch_tidy.csv")
SN<-a[,158:209]
b<-apply(SN,1,function(x){
  ind<-sum(x<1,na.rm=TRUE)/length(x)<0.3
  return(ind)
})
a_0<-a[b,]
write.csv(a_0,"G:/test/output.csv",row.names = F)

## Grades filter(grades has: 'A' 'B' 'C' 'D' '')
a<-read.csv("G:/test/output.csv")
grades<-a[,106:157]
b<-apply(grades,1,function(x){
  ind<-sum(x=='')/length(x)<=0.9
  return(ind)
})
a_0<-a[b,]
write.csv(a_0,"G:/test/output1.csv",row.names = F)

## Height filter(>1000)
a <- read.csv("~/temp/output2_svf_0.1valid.csv", stringsAsFactors = F)
height <- a[, 22:41]
b<-apply(height,1,function(x){
  ind<-mean(x, na.rm=TRUE)/length(x)>1000
  return(ind)
})
a_0<-a[b,]

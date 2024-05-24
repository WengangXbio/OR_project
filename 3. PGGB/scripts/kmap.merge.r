df=read.table("kmap1",head=F)
df1=df[order(df$V2,df$V3),]
fna=df1[1,1]
outdf=data.frame(fna=0,seq=0,start=0,end=0,chr=0)
for(seq in unique(df1$V2)){
chr=as.character(df1[which(df1$V2==seq),5][1])
start=0
end=df1[which(df1$V2==seq),4][1]
tdf=df1[which(df1$V2==seq),]
for(i in 1:(nrow(tdf)-1)){
t1=tdf[i,5]
t2=tdf[i+1,5]
if(t1==t2){
end[length(end)]=c(tdf[i+1,4])
}else{
start=c(start,tdf[i+1,3])
end=c(end,tdf[i+1,4])
chr=c(chr,as.character(tdf[i+1,5]))
}
}
outdf=rbind(outdf,data.frame(fna,seq,start,end,chr))
}
outdf=outdf[-1,]
write.table(outdf,"kmap1.merge",quote=F,sep = "\t",col.names=F, row.names=F)

library(dplyr)
library(gtools)
library(readr)

# from https://gist.github.com/shujishigenobu/1858458 by shujishigenobu
N50 <- function(x) {
    x.sorted <- rev(sort(x))
    return(x.sorted[cumsum(x.sorted) >= sum(x.sorted)*0.5][1])
}

args=commandArgs(T)
statsfn = args[1]
tsv_path = args[2]

st=read.table(statsfn, sep="\t", head=T, row.names=NULL, na=c("-"))
st$alen <- as.numeric(st$alen) # wtf is R

samplename <- st[1,]$samplename
a=st %>% group_by(refname) %>% 
         summarise(n=n(),
                   bases=sum(alen),
                   meanlen=mean(alen),
                   medianlen=median(alen),
                   maxlen=max(alen),
                   N50=N50(alen))

sources=read.table("../ref.cfg", comment="#", sep="\t", head=T, row.names=NULL, na=c("-"))
b=inner_join(a, sources, by="refname")

total=sum(b$bases)
c=b %>% mutate(Cov = bases/(size*1e6)) %>%
        mutate(Perc = bases/total * 100) %>%
        mutate(FoldChange = foldchange(Perc, abundance)) %>%
        mutate(samplename=samplename) %>%
        arrange(Perc)
   
write_tsv(c %>% select(samplename, refname, n, bases, meanlen, medianlen, N50, maxlen, Cov, Perc, abundance, FoldChange, celltype), tsv_path)

#print(c%>%select(refname,bases,N50,Cov,Perc,abundance,FoldChange,celltype))

#p=ggplot(c, aes(x=Genome, y=FoldChange, fill=CellType)) + geom_bar(stat="identity")
#ggsave(paste0(prefix,"_","foldchange.pdf"), p)


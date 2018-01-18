# Calculated FAIME profile of a single sample from its annotation
# Details reference of this method can be found from an article:
# http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002350
# r, a named gene annotation list,
# s, a named rank number of a gene expression order list from  gene expression profile in a single samples

## detailed data structure

# r data structure : annotation list
  #GO1: gene1,gene2,gene3...
  #GO2: gene4,gene5,gene7...
  #...
  #GOn: geneN,geneM,geneJ...

#s is vector with expression rank number, and named by gene name
  #gene1:a,gene2:b,geng3:c...geneN,x
  #like, s=c(1:3);names(s)=c("a","b","c"), abc are gene names, 123 are expression rank order

#output Fscore : the same structure as s

run.FAIME<-function(r,s){
  # get w score  which means the weighted expression values descripted in the text
  # s is the total rank number mentioned above
  wscore<-function(gene,s){
    w<-s[gene]*exp(-s[gene]/length(s))
    return(w)
  }
  #get all wscore of genes in a single sample
  wlist<-c()
  for(i in 1:length(s)){
    wlist<-c(wlist,wscore(names(s)[i],s))
  }
  #asign names from s
  names(w)<-names(s)

  #function
  #w generated from above expression
  get.Fsocre.of.individual.geneset<-funtion(gene.single.set,w){
    NC.GO<-mean(w[names(w)%in%gene.single.set])
    NC.G.GO<-mean(w[!names(w)%in%gene.single.set])
    return(NC.GO-NC.G.GO)
  }
  #gene set names
  c<-names(r)
  Fscore<-c()

  #get Fscore list
  for (i in 1:length(c)){
    gene.single.set<-r[[i]]
    Fscore<-c(Fscore,get.Fsocre.of.individual.geneset(gene.single.set,w))
  }
  #add names
  names(Fscore)<-c



  return(Fscore)

}

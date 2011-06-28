library(MASS)
library(seqinr)


slidingwindowplot <- function(windowsize, inputseq)
{
   print(inputseq)
   starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
   n <- length(starts)    # Find the length of the vector "starts"
   GCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
   for (i in 1:n) {
        chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
        GC <- GC(chunk)
        print(GC)
        GCs[i] <- GC
   }
   plot(starts,GCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}


slidingwindowfreq <- function(windowsize, inputdir, binsize)
{
   GCs <- numeric()
   v <- 0 
   files <- list.files(inputdir,full.names=TRUE)

   for (f in files) {
       singleseqvec <- read.fasta(file = f)
       singleseqnuc<-singleseqvec[[1]]
       lengthseq = length(singleseqnuc) 

       if (lengthseq > windowsize) {

            starts <- seq(1, length(singleseqnuc)-windowsize, by = windowsize)
            n <- length(starts)    # Find the length of the vector "starts"
     
            for (i in 1:n) {
                 chunk <- singleseqnuc[starts[i]:(starts[i]+windowsize-1)]
                 GC <- GC(chunk)
                 GCs[v] <- GC
                 v <- v+1
            }
       }
   } 

   bins <- seq(0,1,by=binsize)
   hist1<-hist(GCs,breaks=bins,plot=FALSE, xlab='GC content %', ylab='% of 1kb bins with this GC content' )
   hist1$counts <- hist1$counts*100/sum(hist1$counts)
   plot(hist1)
}



gcvssize <- function(inputdir)
{
   GCs <- numeric()
   lengths <- numeric()
   v <- 0 
   files <- list.files(inputdir,full.names=TRUE)

   for (f in files) {
       singleseqvec <- read.fasta(file = f)
       singleseqnuc<-singleseqvec[[1]]
       lengthseq = length(singleseqnuc) 

       GC <- GC(singleseqnuc)
       GCs[v] <- GC
       lengths[v] <- lengthseq
       v <- v+1

   } 

   plot(lengths, GCs, ylab="GC content", xlab="Contig length")
}


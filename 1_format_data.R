##
# Functions for formatting source database files
# 2017-10-14
# Author: Song
##

Cut_4DGdata = function()
{
  setwd(YOUR-WORKING-DIRECTORY)
  d<-readLines("4DGenome_HomoSapiens_hg19.txt")
  #d<-read.table('1IMR90_all.txt',sep = '\t')
  d<-unlist(strsplit(d, "\t", fixed = TRUE))
  e<- matrix(d,nrow=3095882,ncol=15,byrow = T)
  
  a<-c(1,2,3,4,5,6,10,11,12,15)
  f<-e[,a]
 
  celltype<-c('IMR90','A549','K562','MCF7','HELA','IPS6.9','IPS19.11','H1ESC')
  
  FILTERBYcelltype<-f[which(f[,7] %in% celltype),]
  
  FILTERBYcelltype<- FILTERBYcelltype[ order(as.integer( FILTERBYcelltype[,2])), ]
  write.table (FILTERBYcelltype, 'All_FILTERBYcelltype_hg19.txt', sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  
  chrtype<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
  
  f<-FILTERBYcelltype[which(FILTERBYcelltype[,1] %in% chrtype),]
  f<- f[ order( f[,1],as.integer( f[,2]),as.integer( f[,3])), ]
  f<-unique(f)
  write.table (f, 'chr11_4DG-hg19.txt', sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  
  C1<-f[,1:3]
  write.table (C1, 'A_chr11_4DG-hg19.txt', sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  #then go to liftOver and generate A_chr11_4DG-hg38.bed!
  
  err<-read.table('A_chr11.err',sep ="\t")
  err <-as.matrix(err)
  #err <-unique(err)
  C1_hg38<-read.table('A_chr11_4DG-hg38.bed',sep ="\t")
  C1_hg38 <-as.matrix(C1_hg38)
  
  j<-1
  for(i in 1:720161)
  {
    if(f[i,1]==err[j,1] & as.integer( f[i,2])==as.integer(err[j,2]) & as.integer(f[i,3])==as.integer(err[j,3]) )
      {
        f[i,1] =0
        f[i,2]=0
        j<-j+1
    }
  }
  e<-f[which(f[,1] != '0'),]
  e[,1:3]<-C1_hg38
  
  C2<-e[,4:6]
  write.table (C2, 'B_chr11_4DG-hg19.txt', sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  
  
  err<-read.table('B_chr11.err',sep ="\t")
  err <-as.matrix(err)
  #err <-unique(err)
  C1_hg38<-read.table('B_chr11_4DG-hg38.bed',sep ="\t")
  C1_hg38 <-as.matrix(C1_hg38)
  j<-1
  for(i in 1:719955)
  {
    if(e[i,4]==err[j,1] & as.integer( e[i,5])==as.integer(err[j,2]) & as.integer(e[i,6])==as.integer(err[j,3]) )
    {
      e[i,1] =0
      e[i,2]=0
      j<-j+1
    }
  }
  f<-e[which(e[,1] != '0'),]
  f[,4:6]<-C1_hg38
  write.table (f, 'chr11_4DG-hg38.txt', sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)

}

Create_TFdata = function()
{
  setwd(YOUR-WORKING-DIRECTORY)
  TF<-read.table('1-TF_human_data_information.txt',sep = '\t')
  TF<-as.matrix(TF)
  a<-c(2,4,7,8)
  b<-c(1:6,9)
  A<-TF[,a]
  
  celltype<-c('IMR90','A549','K562','MCF-7','HeLa','H1','IPS6.9','IPS19.11')
  
  TFbyCT<-A[which(A[,2] %in% celltype),]
 
  for(j in 1:8)
  {
     A549<-TFbyCT[which(TFbyCT[,2] == celltype[j]),]
    flag<-1
  for(i in 1:nrow(A549))
  { 
    e<-read.table(A549[i,4])
    e<-as.matrix(e)
    e[,4]<-A549[i,2]
    e[,5]<-A549[i,1]
    e[,6]<-A549[i,3]
    
    e<-e[,b]
    if(flag==1)
       { ALL<-e
         flag<-0}
      else
      ALL<-rbind(ALL,e)
  }
    for(k in 1:22)
    {
    x <- paste("chr",k,sep = "")
    
     write.table (ALL[which(ALL[,1] == x),], paste(x,'H1_all.txt',sep = ""), sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
    }
  }
     
  ALL<-read.table('H1_all.txt',sep = '\t')
  a<-ALL[,6]
  a<-as.matrix(a)
  a<-unique(a)
  e<-c(a,e)
  e<-unique(e)
  
 
  chrtype<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
  setwd("E:/Database/new_CSDB/A549")
  ALL<-read.table('A549All.txt',sep = '\t')
  ALL<-as.matrix(ALL)
  #ALL<- ALL[ order(ALL[,1],as.integer( ALL[,2])), ]
  #ALL<-ALL[which(ALL[,1] %in% chrtype),]
  A<-ALL
  b<-ALL[,1:2]
  b[,1]<-as.integer(as.integer( ALL[,2])/2000)
  b[,2]<-as.integer(as.integer( ALL[,3])/2000)
  ALL[,2:3]<-b
  
  ALL<- ALL[ order(ALL[,1],as.integer( ALL[,2]),ALL[,6]), ]
  for(i in nrow(ALL):2)
  {
    if(ALL[i,6]==ALL[i-1,6] & ALL[i,1]==ALL[i-1,1] & ALL[i,2]==ALL[i-1,2] &ALL[i,3]==ALL[i-1,3])
    {
      ALL[i,2]<-0
    }
  }
  ALL<-ALL[which(ALL[,2] != 0),]
  for(k in 1:22)
  {
    x <- paste("chr",k,sep = "")
    
    write.table (ALL[which(ALL[,1] == x),], paste(x,'A549All.txt',sep = ""), sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  }
  write.table (ALL[which(ALL[,1] == 'chrX'),], paste('chrX','A549All.txt',sep = ""), sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  
  write.table (ALL[which(ALL[,1] == 'chrY'),], paste('chrY','A549All.txt',sep = ""), sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
   }

Combin_TFRange= function()
{
  library(stringr)
  setwd(YOUR-WORKING-DIRECTORY)
  ALL<-read.table('K562All.txt',sep = '\t')
  
  ALL<-as.matrix(ALL)
  for(n in 1:22) 
  {
  n<-'Y'
  
  A549<-read.table(paste('chr',n,'MCF-7All.txt',sep = ""),sep = '\t')
  A549<-as.matrix(A549)
  #A549<-ALL[which(ALL[,1]==paste('chr',n,sep = "")),]
  A549<- A549[ order(A549[,6],A549[,2]), ]
  A549[,1]<-floor(as.integer( A549[,2])/2000+1)
  A549[,4]<-floor(as.integer( A549[,3])/2000+1)
  A549<-cbind(A549,A549[,5])
  A549[,8]<-''
  i=1
  for(j in 2:nrow(A549))
  {
    if(A549[i,6] == A549[j,6] & A549[i,1] == A549[j,1] & A549[i,4] == A549[j,4])
    { 
        if(A549[i,8]=='')
        {
          A549[i,8]<-paste(str_trim(A549[j,2]),str_trim(A549[j,3]),A549[j,5],paste('Confidence:Qvalue:',str_trim(A549[j,7]),sep = ""),sep = ";")
          A549[j,]<-''
        }
        
      else
      {
        A549[i,8]<-paste(A549[i,8],paste(str_trim(A549[j,2]),str_trim(A549[j,3]),A549[j,5],sep = ";"),sep = "|")
        A549[j,]<-''
      }
    }
    else
      i=j
  }
  A549[,1]<-paste('chr',n,sep = "")
  A549[,4]<-'MCF7'
  write.table (A549[which(A549[,2]!=''),], file =paste(n,'MCF7.txt',sep = ""), sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  }
  
  for(n in 1:22) 
  {
    write.table (ALL[which(ALL[,1]==paste('chr',n,sep = "")),], file =paste(n,'-ADGenome.txt',sep = ""), sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  }  
}

Create_bins= function()
{
  library(stringr)
  for(k in 1:22)
  {
    g<-read.table(paste(k,'-gene-hg38.txt',sep=''),sep='\t')
    g<-as.matrix(g)
    
    a<-nrow(g)
    b<-floor(as.integer( g[a,5])/200000+1)
    f1<-matrix(paste('chr','Y',sep=''),b,5)
    for(i in 1:nrow(f1))
    {
      f1[i,2]<-paste('Bin',i,sep = "")
      f1[i,3]<-(i-1)*200000+1
      f1[i,4]<-i*200000
    }
    f1[,5]<-''
    D<-f1[1,]
    D<-as.matrix(D)
    for(i in 1:nrow(g))
    {
      for(j in floor(as.integer( g[i,4])/200000+1):floor(as.integer( g[i,5])/200000+1))
      {
        C<-f1[which(f1[,2]==paste('Bin',j,sep = "")),]
        C<-as.matrix(C)
        C[5]<-paste(str_trim(g[i,4]),str_trim(g[i,5]),str_trim(g[i,9]),sep=';')
        D<-cbind(D,C)
      }
    }
    
    D<-t(D)   #转置
    A<-rbind(f1,D)
    A<- A[ order(as.integer(A[,3])), ]
    j=1
    for(i in 2:nrow(A))
    {
      if( A[i,2]==A[j,2] )
      {
        if(A[j,5]=='')
          A[j,5]<-A[i,5]
        else
          A[j,5]<-paste(A[j,5],A[i,5],sep = "|")
      }
      else
        j=i
    }
    j=1
    for(i in 2:nrow(A))
    {
      if( A[i,2]==A[j,2] )
      {
        A[i,5]<-A[j,5]
      }
      else
        j=i
    }
    A<-unique(A)
    write.csv(A, paste('h','Y','-ADGenome-hg38.csv',sep=''))
  }
}

Create_Range= function()
{
  h<-0
  for(i in 1:22)
  {
    g<-read.table(paste(i,'-ADGenome-h38.txt',sep = ""),sep='\t')
    g<-as.matrix(g)
    a<-c(1,2,3)
    g1<-g[,a]
    h<-rbind(h,g1)
    b<-c(4,5,6)
    g2<-g[,b]
    h<-rbind(h,g2)
    h<-unique(h)                     #把重复的元素去掉
  }
  
  g<-cbind(h[,1],str_trim(h[,2]),str_trim(h[,3]))
  for(i in 1:22)
  {
    write.csv(g[which(g[,1]==paste('chr',i,sep = "")),], file =paste(i,'gene-hg38.csv',sep = ""))
  }
  
  for(i in 1:22)
  {
    write.table (g[which(g[,1]==paste('chr',i,sep = "")),], file =paste(i,'-gene-hg38.txt',sep = ""), sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  }
  for(i in 1:22)
  {
    g<-read.csv(paste(i,'gene-hg38.csv',sep = ""))
    g<-as.matrix(g)
    g<-unique(g)
    write.csv(g, file =paste(i,'-gene-hg38.csv',sep = ""))
  }
  write.csv(h, 'gene-hg38.csv')
}



Create_iRefInedx= function()
{
  setwd(YOUR-WORKING-DIRECTORY)
  f<-read.csv('TF-TF-2.csv',sep=',')
  f<-as.matrix(f)
  TF<-read.csv('TF-new.csv',sep=',')
  TF<-as.matrix(TF)
  
  a<-c(1,2,5,6,7,9,15,51,53)
  f<-f[,a] 
  y<-f[1:2,]
  {
  for(i in 1:261)
  {
    #x<-f[which(regexpr(TF[i],f[,1])[1]>1),]
    x<-f[1:2,]
    for(j in 1:nrow(f))
    {
      if(regexpr(TF[i],f[j,3])[1]>1)
      {
        x<-rbind(x,f[j,])
      }
    }
    if(nrow(x)>2)
    {
      x[,3]=TF[i]
      x<-x[-1,]
      x<-x[-1,]
      y<-rbind(y,x)
   }
  }
  
  y<-y[-1,]
  y<-y[-1,]
  f<-y[1:2,]
  for(i in 1:261)
  {
    #x<-f[which(regexpr(TF[i],f[,1])[1]>1),]
    x<-y[1:2,]
    for(j in 1:nrow(y))
    {
      if(regexpr(TF[i],y[j,4])[1]>1)
      {
        x<-rbind(x,y[j,])
      }
    }
    x[,4]=TF[i]
    if(nrow(x)>2)
    {
     x<-x[-1,]
     x<-x[-1,]
     f<-rbind(f,x)
    }
  }
  f<-f[-1,]
  f<-f[-1,]
  f<-unique(f)
  }
  
  f<-read.csv('TF-TF-comb.csv',sep=',')
  f<-as.matrix(f)
  j<-1
  for(i in 2:8498)
  {
    if(f[i,1]==f[j,1] & f[i,2]==f[j,2])
    {
      f[j,5]=paste( f[j,5],f[i,5],sep = '|')
      f[i,5]<-0
      
      if(f[i,6] %in% f[j,6])
        f[i,6] <-0
      else
        f[j,6]=paste( f[j,6],f[i,6],sep = '|')
    }
    else
      j<-i
  }
  
  j<-1930
  k<-1
  for(i in 1931:4911)
  {
    if(f[i,1]==f[j,1] & f[i,2]==f[j,2])
    {
      if(k==1)
         {
          k<-k+1
          f[j,3]=paste(paste( "Rel1:",f[j,3],sep = ''),paste( "Rel",k,":",f[j,3],sep = ''),sep = '|')
        }
      else
      {
        k<-k+1
        f[j,3]=paste(f[j,3],paste( "Rel",k,":",f[j,3],sep = ''),sep = '|')
      }
      f[i,3]<-0
    }
    else
    {
      k<-1
      j<-i
    }
      
  }
  
  f<-f[which(f[,3] != 0),]
  
  write.csv(f, 'TF-TF-2.csv')
}

Cut_LncRNAdata = function()
{
  d<-readLines("high_hg38_uniq.txt")
  d<-unlist(strsplit(d, "\t", fixed = TRUE))
  e<- matrix(d,nrow=270477,ncol=24,byrow = T)   #nrow=270477
  
  e<-as.matrix(d)
  a<-c(1,2,8,9,14,18,20,22,23,24)
  f<-e[,a]
  name<-f[,6]
  name <- unique(name)
  celltype<-c('H1')
  #'MCF7 cell','HeLa'
  f<-f[which(f[,6] == 'H1'),]
  f[,6]<-'H1ESC'
  d<-e[,a]
  b<-d[which(d[,6] == 'MCF7 cell'),]
  b[,6]<-'MCF7'
  f<-rbind(f,b)
  
  f<- f[ order(as.integer( f[,9])), ]
  write.table (f,'LR_high_hg38-all.txt', sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
  
}

Unique_Range= function()
{
  library(stringr)
  for(i in 1:22)
  {
    f<-read.csv(paste('Y','gene-hg38.csv',sep = ""),sep=',')
    f<-as.matrix(f)
    H<-read.csv(paste('Y','A-Range-hg38.csv',sep = ""),sep=',')
    H<-as.matrix(H)
    H<-H[,-1]
    H<-H[,-1]
    
    H[,1]<-str_trim(H[,1])
    H[,2]<-str_trim(H[,2])
    h<-rbind(f,H)
    h<-unique(h)                     #把重复的元素去掉
    
    
    write.csv(h, file =paste('Y','-gene-hg38.csv',sep = ""))
    
  }
  
}

Unique_TF= function()
{
  TF <- read.csv('TF-TF-2.csv')
  TF <- as.matrix(TF)
  j=1   #把相同片段的 细胞类型值设为相同
  for(i in 2:nrow(f))
  {
    if( f[i,3]==f[j,3] & f[i,2]==f[j,2] & f[i,5]==f[j,5] )
    {
      f[j,8]<- paste(f[j,8],f[i,8],sep='|')
    }
    else
      j=i
  }
  
  j=1   #把相同片段的 细胞类型值设为相同
  for(i in 2:nrow(f))
  {
    if( f[i,3]==f[j,3] & f[i,2]==f[j,2] & f[i,5]==f[j,5] )
    {
      f[i,8]<- f[j,8]
    }
    else
      j=i
  }
  f<-unique(f)
  
  write.csv(f, 'TF-TF-2.csv')
  j=1   #把相同片段的 细胞类型合并
  k<-1;
  for(i in 2:nrow(TF))
  {
    if( TF[i,1]==TF[j,1] & TF[i,2]==TF[j,2] )
    {
      if(k==1)
      { 
        TF[j,3]<- paste (paste("Rel1",TF[j,3],sep=":"), paste("Rel2",TF[i,3],sep=":"),sep="|")
         k<-k+1
      }
      else
      { 
        TF[j,3]<- paste (TF[j,3], paste(paste("Rel",k,sep=""),TF[i,3],sep=":"),sep="|")
        k<-k+1
      }
    }
    else
    {
      j=i
      k=1
     }
  }
  
  j=1   #把相同片段的 细胞类型值设为相同
  for(i in 2:nrow(TF))
  {
    if( TF[i,1]==TF[j,1] & TF[i,2]==TF[j,2]  )
    {
      TF[i,3]<- TF[j,3]
    }
    else
      j=i
  }
  
  TF<-unique(TF)                     #把重复的元素去掉
  
  
  TF<-unique(TF)  
  
  TF<- TF[ order(TF[,3],TF[,4]), ]      #排序
  
  write.csv(TF, 'TF-unique.csv')
}

Unique_TF_Bind= function()
{
  library(stringr)
  A<-read.table('1IMR90_all.txt',sep = '\t')
  A<-as.matrix(A)
  for(n in 1:22)
  {
    #ALL<- A[which(A[,1] == paste('chr','n',sep = "")),]
    ALL<-read.table(paste(n,'MCF7.txt',sep = ""),sep = '\t')
    ALL<-as.matrix(ALL)
    ALL<- ALL[order(ALL[,6],ALL[,2]), ] 
    ALL[,1] = floor(as.integer( ALL[,2])/2000+1) 
    ALL[,4] = floor(as.integer( ALL[,3])/2000+1)
    
    #ALL<-cbind(ALL,ALL[,5])
    #ALL[,7]<-''
    i=1
    for(j in 2:nrow(ALL))
    {
      if(ALL[i,6] == ALL[j,6] & ALL[i,1] == ALL[j,1] & ALL[i,4] == ALL[j,4])
      { 
        if(ALL[i,7]=='')
        {
          ALL[i,7]<-paste(str_trim(ALL[j,2]),str_trim(ALL[j,3]),ALL[j,5],sep = ";")
          ALL[j,]<-''
        }
        
        else
        {
          ALL[i,7]<-paste(ALL[i,7],paste(str_trim(ALL[j,2]),str_trim(ALL[j,3]),ALL[j,5],sep = ";"),sep = "|")
          ALL[j,]<-''
        }
      }
      else
        i=j
    }
    ALL[,1]<-paste('chr',n,sep = "")
    ALL[,4]<-'MCF7'
    
    write.table (ALL[which(ALL[,2]!=''),], file =paste(n,'A-MCF7.txt',sep = ""), sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
    
  }
  
}

Combin_All_CEll_TF_Bind= function()
{
  library(stringr)
  for(n in 1:22)
  {
    n<-1
    A<-read.table(paste(n,'A-A549.txt',sep = ""),sep = '\t')
    A<-as.matrix(A)
    B<-read.table(paste(n,'A-H1.txt',sep = ""),sep = '\t')
    B<-as.matrix(B)
    C<-read.table(paste(n,'A-Hela.txt',sep = ""),sep = '\t')
    C<-as.matrix(C)
    D<-read.table(paste(n,'A-IMR90.txt',sep = ""),sep = '\t')
    D<-as.matrix(D)
    E<-read.table(paste(n,'A-K562.txt',sep = ""),sep = '\t')
    E<-as.matrix(E)
    F<-read.table(paste(n,'A-MCF7.txt',sep = ""),sep = '\t')
    F<-as.matrix(F)
    
    ALL<-rbind(A,B)
    ALL<-rbind(ALL,C)
    ALL<-rbind(ALL,D)
    ALL<-rbind(ALL,E)
    ALL<-rbind(ALL,F)
    
    ALL<-cbind(ALL[,1],ALL)
    ALL<- ALL[order(ALL[,7],ALL[,3]), ]
    
    ALL[,1] = floor(as.integer( ALL[,3])/2000+1) 
    ALL[,2] = floor(as.integer( ALL[,4])/2000+1)
    
    i=1
    for(j in 2:nrow(ALL))
    {
      if(ALL[i,7] == ALL[j,7] & ALL[i,1] == ALL[j,1] & ALL[i,2] == ALL[j,2])
      { 
        if(ALL[i,8]=='')
        {
          ALL[i,8]<-paste(str_trim(ALL[j,3]),str_trim(ALL[j,4]),ALL[j,6],sep = ";")
          ALL[j,]<-''
        }
        
        else
        {
          ALL[i,8]<-paste(ALL[i,8],paste(str_trim(ALL[j,3]),str_trim(ALL[j,4]),ALL[j,6],sep = ";"),sep = "|")
          ALL[j,]<-''
        }
      }
      else
        i=j
    }
    ALL[,1]<-paste('chr',n,sep = "")
    
    write.table (ALL[which(ALL[,2]!=''),], file =paste(n,'aaA-Hela.txt',sep = ""), sep ="\t", row.names =FALSE, col.names =FALSE, quote =FALSE)
    
  }
}

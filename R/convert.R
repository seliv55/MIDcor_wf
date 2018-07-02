
run_convert<-function(infile="midcorout.csv",outfile='smprow'){
# converts Metabolights exchange format into that convenient for visual check
  write("",outfile);
  rada<-read.table(infile, sep=",");   # read experimental data
  tit<-data.frame(lapply(rada[1,], as.character), stringsAsFactors=FALSE)
  rada<- droplevels(rada[-1,])
  print(tit)
        abund<- grep("abundance", tit)[2]  # calculated fractions of isotopologs
        colmet<- grep("Metab", tit)  # column of metabolite name
        colfrg<- grep("atomic", tit)  # carbon positions in the fragment
        coltrac<- grep("tracer", tit)  # carbon positions in the fragment
        colab<- grep("labell", tit)  # carbon positions in the fragment
        colfr<- grep("abundance", tit)[1]  # carbon positions in the fragment
        colcel<- grep("cells", tit)  # carbon positions in the fragment
        coltim<- grep("time", tit)  # carbon positions in the fragment
    levtrac<-levels(rada[,coltrac][drop=T])
    levlpos<-levels(rada[,colab][drop=T])
    levabun<-levels(rada[,colfr][drop=T])
   for(trlev in levtrac){ iso<- 0; abun<- 0 # select tracer
        trc<- subset(rada,trlev==rada[,coltrac])
        if(nchar(trlev)>3){
         llev <- as.character(trc[1,colab])
          a<- strsplit(llev,',')[[1]]
          ara<- length(a)-which(grepl('1',a))
          for(il in 1:length(ara)) iso<- iso+2**ara[il]
          abun <- as.numeric(as.character(trc[1,colfr]))/100.
        }
        write(paste("tracer:",trlev,iso,abun),outfile,append=T)
        levcel<-levels(trc[,colcel][drop=T])
   for(celev in levcel){ # select cells
        cll<- subset(trc,celev==trc[,colcel])
        write(paste("CELLS:",celev),outfile,append=T);
        levmet<- levels(cll[,colmet][drop=T])
   for(mlev in levmet){ # select metabolite
        met<- subset(cll,mlev==cll[,colmet])
        levfr<- levels(met[,colfrg][drop=T]) 
    for(frlev in levfr){ # select fragment
        metfr<- subset(met,frlev==met[,colfrg])
       pCf<-gregexpr("C",frlev)[[1]]; # C-position in formula
       iCb<- as.numeric(substr(frlev,pCf[1]+1,pCf[1]+1))-1
       iCe<- as.numeric(substr(frlev,pCf[2]+1,nchar(frlev)))-1
       nCfrg<- iCe-iCb+1
        write(paste(paste("name:",mlev),iCb,iCe,sep=','),outfile,append=T);
        levtim<- levels(metfr[,coltim][drop=T])
    for(tilev in levtim){ # select fragment
        tot<-data.frame();
        metfrt<- subset(metfr,tilev==metfr[,coltim])
        if(nchar(tilev)>0) write(paste("t=",tilev),outfile,append=T) else
         write("t= 0",outfile,append=T)
        levsmp<- levels(metfrt[,1][drop=T])
     for(smplev in levsmp){ # select sample
        df<- data.frame()
        smp<- subset(metfrt,(smplev==metfrt[,1])&(!grepl('-1',metfrt[,abund-1])))
        dis<- as.numeric(as.character(smp[,abund]))
        df[1,1]<- smp[1,1];
        for(idis in 1:(nCfrg+1)) df[1,1+idis]<- round(dis[idis],4)
    tot<- rbind(tot,df)
      	   }
   write.table(tot,outfile,sep=",",append=T,col.names=FALSE, row.names = F);
          }
         }
        }
       }
      }
}
#   write.table(tot,outfile,sep=",",append=TRUE,col.names=FALSE, row.names = F);

msdlist<-function(trati) { nln<-length(trati)
        ntrati<-lapply(strsplit(trati, " "),as.numeric)
        lnumv<-!lapply(ntrati,is.na)[[1]]
        mdis<-matrix(ncol=length(ntrati[[1]][lnumv]),nrow=nln,0)
        mdis[1,]<-ntrati[[1]][lnumv]
      if(nln>1)for(i in 2:nln){
        lnumv<-!lapply(ntrati,is.na)[[i]]
        mdis[i,]<-ntrati[[i]][lnumv]
      } 
        mval<-round(apply(mdis,2,mean),4)
        sdv<-round(apply(mdis,2,sd))
 return(list(mval,sdv))}
 
 getim<-function(star="cntn",pat="t= "){
  roti<-grep(pat,cntn)
  a<-star[roti]
  aa<-sub(pat,'',a)
  aaa<-'0';
  for(i in aa) if(!(i %in% aaa)) aaa<-c(aaa,i)
  return()
 }
 
 isoform<-function(isofi='smprow'){
        cntn<-readLines(isofi)
        rowti<-grep('t= ',cntn)
        a<- sub('t= ','',cntn[rowti])
        b<-paste(sort.int(unique(c('0',a))),collapse=' ')
        inti<- paste('time(h):',b,'-1',collapse=' ')
        rowtrac<-grep('trac',cntn)  #select samples of the same tracer
        rowtrac<-c(rowtrac,length(cntn)+1)
    for(itr in 1:(length(rowtrac)-1))
      if(nchar(cntn[rowtrac[itr]])>15){
       atr<-cntn[rowtrac[itr]:(rowtrac[itr+1]-1)]
       trspl<- strsplit(atr[1],' ')[[1]][2]
       trspl<-tail(strsplit(trspl,'-')[[1]],n=1)
       rowcel<- grep('CELL',atr)  #select samples of the same cell type
       rowcel<- c(rowcel,length(atr)+1)
      for(icel in 1:(length(rowcel)-1)){
        acel<- atr[rowcel[icel]:(rowcel[icel+1]-1)]
        celspl<- strsplit(acel[1],' ')[[1]]
         fi<-paste(celspl[2],'-',trspl,sep='')
          write(inti,fi);
         rowname<- grep('name',acel)  #select metabolite
         rowname<- c(rowname,length(acel)+1)
       for(ina in 1:(length(rowname)-1)){
          aname<- acel[rowname[ina]:(rowname[ina+1]-1)]
          naspl<- strsplit(aname[1],' ')[[1]]
          write(aname[1],fi,append=T);
          rowti<- grep('t=',aname)  #select time
          rowti<- c(rowti,length(aname)+1)
         for(iti in 1:(length(rowti)-1)){
           ati<- aname[rowti[iti]:(rowti[iti+1]-1)]
          write(ati[1],fi,append=T);
           tispl<- strsplit(ati,',')
           dis<- matrix(nrow=length(ati)-1,ncol=length(tispl[[2]])-1)
           for(i in 1:nrow(dis)) dis[i,1:ncol(dis)]<-0.01*as.numeric(tispl[[i+1]][2:(ncol(dis)+1)])
           mdis<- apply(dis,2,mean)
           sdis<- apply(dis,2,sd)
          write("mean:",fi,append=T);
          write.table(round(t(mdis),3),fi,append=T, col.names = F, row.names = F);
          write("sd:",fi,append=T);
          write.table(round(t(sdis),3),fi,append=T, col.names = F, row.names = F);
           }
         }
         write(atr[1],fi,append=T);
       }
     }
        
}

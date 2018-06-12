
run_convert<-function(infile="midcorout.csv",outfile='smprow'){
  write("",outfile);
  rada<-read.table(infile, sep=",");   # read experimental data
  tit<-data.frame(lapply(rada[1,], as.character), stringsAsFactors=FALSE)
  rada<- droplevels(rada[-1,])
  print(tit)
        abund<- grep("abundance", tit)[2]  # calculated fractions of isotopologs
        colmet<- grep("Metab", tit)  # column of metabolite name
        colfrg<- grep("atomic", tit)  # carbon positions in the fragment
        coltrac<- grep("tracer", tit)  # carbon positions in the fragment
        colcel<- grep("cells", tit)  # carbon positions in the fragment
        coltim<- grep("time", tit)  # carbon positions in the fragment
    levtrac<-levels(rada[,coltrac])
   for(trlev in levtrac){ # select tracer
        trc<- subset(rada,trlev==rada[,coltrac])
        if(nchar(trlev)>3) write(paste("TRACER:",trlev),outfile,append=T)
          else write("TRACER: N/A",outfile,append=T)
        levcel<-levels(trc[,colcel][drop=T])
   for(celev in levcel){ # select cells
        cll<- subset(trc,celev==trc[,colcel])
        write(paste("CELLS:",celev),outfile,append=T);
        levmet<- levels(cll[,colmet][drop=T])
   for(mlev in levmet){ # select metabolite
        met<- subset(cll,mlev==cll[,colmet])
        levfr<- levels(met[,colfrg][drop=T]) 
    for(frlev in levfr){ # select fragment
        tot<-data.frame();
        metfr<- subset(met,frlev==met[,colfrg])
        write(paste("Metabolite:",mlev,"Frg:",frlev),outfile,append=T);
        levtim<- levels(metfr[,coltim][drop=T])
    for(tilev in levtim){ # select fragment
        metfrt<- subset(metfr,tilev==metfr[,coltim])
        write(paste("Incubation",tilev),outfile,append=T);
        levsmp<- levels(metfrt[,1][drop=T])
     for(smplev in levsmp){ # select sample
        df<- data.frame()
        smp<- subset(metfrt,smplev==metfrt[,1])
        dis<- as.numeric(as.character(smp[,abund]))
        df[1,1]<- smp[1,1];
        for(idis in 1:nrow(smp)) df[1,1+idis]<- round(dis[idis],4)
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
   


correct<-function(chast,fn1,mdcor){
    md<-substr(mdcor,1,2)
  onemet<-convert(chast) 
  id<-onemet[[1]]; mm<-onemet[[2]] # matrix of signal intensities
   nC<-onemet[[3]];  nfrg<-onemet[[4]]; nSi=onemet[[5]];  nS<-onemet[[6]]
  labmet<-onemet[[7]]
   numc<- ncol(mm);  nmass<- nfrg+1; 
   if(numc==nmass) {mm<-cbind(mm,mm[,numc]); numc<- ncol(mm);}
    nln<-length(id); mdful<-mm

      mmteor<-mtr(nfrg,numc,nC,nSi,nS) # theoretic distribution
         
# mass fractions
   fr<-mdistr(nfrg,mm,mmteor,nln) # write mass fractions without correction:
 write("*** MID for each injection, corrected only for natural 13C, 29,30Si, 33,34S ***",fn1,append=T)
  write.table(cbind(id,round(fr,5)),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
        write("\n",fn1,append=TRUE);

# correction factor
     corr<-numeric(nmass);
     for(i in 1:nln) {if(nchar(labmet[i])<5) { print(paste(id[i], labmet[i]," len=",nchar(labmet[i]))); break;}}
       corr<-mm[i,1:nmass]-mmteor[1,1:nmass]    # correction factor
       
 if(md=="va") { for(ii in 1:9) {
# correction for EACH INJECTION
    tmp<-mdful;
    for(j in 1:nln)  for(k in 1:(nfrg+1))  for(i in 1:(nmass+1-k)) tmp[j,i+k-1]<-tmp[j,i+k-1]-corr[i]*(fr[j,k]);
      fr<-mdistr(nfrg,tmp,mmteor,nln);
        }}
  if(md=="co") {     # ALL INJECTIONS:
                 for(j in 1:nln)  mdful[j,1:nmass]<-mdful[j,1:nmass]-corr;
     fr<-100*round(mdistr(nfrg,mdful,mmteor,nln),5)
               }
write("*** All samples fully corrected **",fn1,append=T)
  write.table(cbind(id,round(fr,2)),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F)
        write("\n",fn1,append=TRUE);
#
     i<-1 # Statistics for SUMMED INJECTIONS
#  while(i<nlnsum){ cnlab<-substr(idsum[i],1,6); k<-i
#    while(substr(idsum[k+1],1,6)==cnlab) k<-k+1
#    # mean for SUMMED INJECTIONS
#      if(k>i) {sredn<-frsum[1,]; sid<-paste("mean",substr(idsum[i],1,6),sep="_")
#                 for(j in 1:nmass) sredn[j]<-mean(frsum[i:k,j])
#                 for(j in 1:nmass) if(sredn[j]<0.) sredn[j]<-0.
#                 srsum<-sum(sredn[1:nmass])
#               sredn[1:nmass]<-sredn[1:nmass]/srsum # renormalization
#    # SD for SUMMED INJECTIONS
#               otkl<-frsum[1,]; oid<-paste("sd",substr(idsum[i],1,6),sep="_");
#                 for(j in 1:nmass) otkl[j]<-sd(frsum[i:k,j]);
#        frsum<-rbind(frsum[1:k,],sredn,otkl,frsum[(1+k):nrow(frsum),]);
#        idsum<-c(idsum[1:k],sid,oid,idsum[(1+k):length(idsum)]);
#    write.table(t(c(sid,nfrg,round(sredn[1:(nfrg+1)],5))),fn2,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
#    write.table(t(c(oid,nfrg,round(otkl[1:(nfrg+1)],5))),fn2,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
## write.table(t(c("*Correction: ",round(corr,5))),fn2,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
#                  i<-k+2;}
#          i<-i+1; nlnsum<-nrow(frsum);}
# write data:
# write("*** Statistics, samples fully corrected **",fn1,append=TRUE)
#    write.table(cbind(idsum,round(frsum,5)),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
 write.table(t(c("*Correction: ",round(corr,5))),file=fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
        write("\n",file=fn1,append=TRUE);
         return(list(id,fr,nfrg))
}

run_midcor<-function(infile="../readCDF/RaMID/ramidout.csv", outfile="midcorout.csv",mode="con"){
  fn1<-paste(infile,"_c",sep="");	
  write("",fn1);
  write("",outfile);
  rada<-read.table(infile, sep=",");   # read experimental data
  tit<-data.frame(lapply(rada[1,], as.character), stringsAsFactors=FALSE)
        abund<- grep("abundance", tit)[2]  # calculated fractions of isotopologs
        colmet<- grep("Metab", tit)  # column of metabolite name
        colfrg<- grep("atomic", tit)  # carbon positions in the fragment
        tot<-data.frame(); #for rada corrected
    Mtb<-levels(rada[,colmet])
    Mtb<-Mtb[!grepl("name",Mtb)][drop=T]
   for(i in 1:length(Mtb)){
        write(paste("----",Mtb[i]," ---"),fn1,append=TRUE)
        met<- subset(rada,Mtb[i]==rada[,colmet]) # select metabolite
        levmet<- levels(met[,colfrg][drop=T]) # fragments of selected met
    for(ilev in 1:length(levmet)){
   metfr<- subset(met,levmet[ilev]==met[,colfrg])
   metfr<-rbind(tit,metfr)
    res<-correct(metfr,fn1,mode)
    id<- res[[1]]
    mcor<- res[[2]]
    for(j in 1:length(id)){
    smetfr<- subset(metfr,metfr[,1]==id[j])
    mcol<- c(0,mcor[j,])
    mcol<- mcol[-length(mcol)]
    racor<- cbind(smetfr[-abund],mcol[1:nrow(smetfr)])
    tot<- rbind(tot,racor)
          }
         }
        }
   write.table(rada[1,],outfile,sep=",",append=F,col.names=FALSE, row.names = F);
   write.table(tot,outfile,sep=",",append=TRUE,col.names=FALSE, row.names = F);
   
        colcel<- grep("cell", tit)  # column of cell type(conditions)
   cells<-levels(rada[,colcel])
   cells<-cells[!grepl("cells",cells)][drop=T]
   
        coltrac<- grep("tracer", tit)  # column of tracer used
   tracer<-levels(rada[,coltrac])
   tracer<-tracer[grepl("C13",tracer)][drop=T]
   wd<-strsplit(outfile,'/')[[1]]
   
   for(i in 1:length(cells)) {
   conds<-subset(tot,(grepl(cells[i],tot[,2])))
   ficond=paste(cells[i])
   if(grepl('/data/',outfile)) ficond<-paste('/data/',ficond,sep="")
   
   write.table(rada[1,],ficond,sep=",",append=F,col.names=FALSE, row.names = F);
   write.table(conds,ficond,sep=",",append=TRUE,col.names=FALSE, row.names = F);
   }
   fiso='smprow'
    run_convert(infile=outfile,outfile=fiso)
    isoform(isofi=fiso)
   return(tot) }


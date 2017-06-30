  
correct<-function(chast,fn1,mdcor){
  md<-substr(mdcor,1,2);
  onemet<-convert(chast,2); #onemet=c("CDF file","spectra all inj","line#","nmet","nC","nfrg","spectra summed inj","id sum")
  id<-onemet[[1]]; mm<-onemet[[2]] # signal intensities for mass range for ALL inj
   iln<-onemet[[3]]; nmet<-onemet[[4]]; nC<-onemet[[5]];
  nfrg<-onemet[[6]]; nSi=onemet[[7]] # signal intensities for mass range for SUMMED inj
   nS<-0;
   numc=ncol(mm);  nmass=nfrg+1; if(numc==nmass) {mm<-cbind(mm,mm[,numc]); numc=ncol(mm);}
    nln<-length(id); mdful<-mm

      mmteor<-mtr(nfrg,nmass,nC,nSi,nS) # theoretic distribution
         
# mass fractions
   fr<-mdistr(nfrg,mm,mmteor,nln) # write mass fractions without correction:
 write("*** MID for each injection, corrected only for natural 13C, 29,30Si, 33,34S ***",fn1,append=T)
  write.table(cbind(id,round(fr,5)),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
        write("\n",fn1,append=TRUE);

# correction factor
     corr<-numeric(nmass);
     for(i in 1:nln) {if(grepl("nat",id[i])||grepl("com",id[i])||grepl("com",id[i])||grepl("COLD",id[i])||grepl("Cold",id[i])) break;}
       corr<-mm[i,1:nmass]-mmteor[1,1:nmass]    # correction factor
       
 if(md=="va") { for(ii in 1:9) {
# correction for EACH INJECTION
    tmp<-mdful;
    for(j in 1:nln)  for(k in 1:(nfrg+1))  for(i in 1:(nmass+1-k)) tmp[j,i+k-1]<-tmp[j,i+k-1]-corr[i]*(fr[j,k]);
      fr<-mdistr(nfrg,tmp,mmteor,nln);
        }}
  if(md=="co") {     # ALL INJECTIONS:
                 for(j in 1:nln)  mdful[j,1:nmass]<-mdful[j,1:nmass]-corr;
     fr<-mdistr(nfrg,mdful,mmteor,nln)
               }
write("*** All samples fully corrected **",fn1,append=T)
  write.table(cbind(id,round(fr,5)),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F)
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
         return(list(fr,nfrg,iln))
}

row2col<-function(res,metss,tot,isoname){
  nstrok<-nrow(res[[1]])
  i<-1; newcol<-numeric()
    for(j in 1:nstrok){
     while(!grepl("C0",metss[i,isoname])) {i=i+1;}
      for(k in 1:(res[[2]]+1)) {newcol[i]<-res[[1]][j,k]*100.; i=i+1; }
              }
           newcol<-round(newcol,3);   if(length(newcol)<nrow(metss)) newcol<-c(newcol,0)
        metss<-cbind(metss[-c(1),],newcol[-c(1)]); tot<-rbind(tot,metss)
    return(tot)}


run_midcor<-function(infile="ramidout.csv", outfile="midcorout.csv",mode="con"){
  fn<-file.path(infile);
  fn1<-paste(fn,"_c",sep="");
   splifi<-strsplit(infile,"/")[[1]]
    if(length(splifi)>1) datadir<-paste("/",splifi[2],"/",sep="") else datadir<-""
  write("",fn1);
  write("",outfile);
  rada<-read.table(fn, sep=",");   # read experimental data
  for(i in 1:ncol(rada)) {
        if(grepl("tracer molecule",rada[1,i])) {coltrac<-i}  # column of tracer used
        if(grepl("cell",rada[1,i])) {colcel<-i}  # column of cell type(conditions)
        if(grepl("time",rada[1,i])) {coltime<-i}  # column of incubation time
         if(grepl("Metab",rada[1,i])) {colmet<-i}  # column of metabolite name
        if(grepl("atomic pos",rada[1,i])) {colfrg<-i}  # carbon positions in the fragment
        if(grepl("signal intensity",rada[1,i])) {inten<-i}  #  column of signal intensity
        if(grepl("isotopol",rada[1,i])) {isoname<-i; } # isotopolog symbols
        if(grepl("abundance",rada[1,i])) {abund<-i}  # calculated fractions of isotopologs
  }
   Mtb<-levels(rada[,colmet])
   Mtb<-Mtb[!grepl("e name",Mtb)][drop=T]
tot<-data.frame(); # group the data for each metabolite and correct for natiral abundance
   for(i in 1:length(Mtb)){
        write(paste("----",Mtb[i]," ---"),fn1,append=TRUE)
   metfr<-rbind(rada[1,],subset(rada,Mtb[i]==rada[,colmet]))
    res<-correct(metfr,fn1,mode)
    tot<-row2col(res,metfr,tot,isoname)
   }
    tot<-tot[-(abund)]
    tot[,inten]<-as.numeric(levels(tot[,inten])[tot[,inten]])
   write.table(rada[1,],outfile,sep=",",append=F,col.names=FALSE, row.names = F);
   write.table(tot,outfile,sep=",",append=TRUE,col.names=FALSE, row.names = F);
   
   cells<-levels(rada[,colcel])
   cells<-cells[!grepl("cells",cells)][drop=T]
   
   tracer<-levels(rada[,coltrac])
   tracer<-tracer[grepl("C13",tracer)][drop=T]
   
   for(i in 1:length(cells)) {
   conds<-subset(tot,(grepl(cells[i],tot[,2])))
   ficond=paste(datadir,cells[i],sep="")
   write.table(rada[1,],ficond,sep=",",append=F,col.names=FALSE, row.names = F);
   write.table(conds,ficond,sep=",",append=TRUE,col.names=FALSE, row.names = F);
   }
    
   return(tot) }


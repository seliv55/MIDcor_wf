  
correct<-function(chast,fn1,fn2,mdcor){
  md<-substr(mdcor,1,2);
  onemet<-convert(chast,2); #onemet=c("CDF file","spectra all inj","line#","nmet","nC","nfrg","spectra summed inj","id sum")
  id=onemet[[1]]; mm=onemet[[2]] # signal intensities for mass range for ALL inj
   iln=onemet[[3]]; nmet<-onemet[[4]]; nC<-onemet[[5]];
  nfrg<-onemet[[6]]; mmsum=onemet[[7]]; idsum=onemet[[8]] # signal intensities for mass range for SUMMED inj
   nSi<-0; nS<-0;
   numc=ncol(mm);  nmass=nfrg+1; nln<-length(id); nlnsum<-length(idsum)

      mmteor<-mtr(nfrg,nmass,nC,nSi,nS) # theoretic distribution
       ef<-sum(mmsum[,1])/sum(mmsum[,2]) # electronic impact factor
         
# normalization for EACH INJECTION
    mdful<-Norm(mm,ef);
# mass fractions
   fr<-mdistr(nfrg,mdful,mmteor,nln) # write mass fractions without correction:
 write("*** MID for each injection, corrected only for natural 13C, 29,30Si, 33,34S ***",fn1,append=T)
  write.table(cbind(id,round(fr,5)),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
        write("\n",fn1,append=TRUE);

# normalization for SUMMED INJECTIONS of each replicate
         gcmsum<-Norm(mmsum,ef);
# mass fractions
      frsum<-mdistr(nfrg,gcmsum,mmteor,nlnsum);
 write("*** Summed injections for each plate, corrected only for natural 13C, 29,30Si, 33,34S **",fn1,append=T)
  write.table(cbind(idsum,round(frsum,5)),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F)
        write("\n",fn1,append=TRUE);
# correction factor
     corr<-numeric(nmass);
     for(i in 1:nlnsum) {if(grepl("nat",idsum[i])||grepl("com",idsum[i])) break;}
       corr<-gcmsum[i,1:nmass]-mmteor[1,1:nmass]    # correction factor
       
 if(md=="va") { for(ii in 1:9) {
# correction for SUMMED INJECTIONS
    tmp<-gcmsum;
    for(j in 1:nlnsum)  for(k in 1:(nfrg+1))  for(i in 1:(nmass+1-k)) tmp[j,i+k-1]<-tmp[j,i+k-1]-corr[i]*(frsum[j,k]);
      frsum<-mdistr(nfrg,tmp,mmteor,nlnsum);
# correction for EACH INJECTION
    tmp<-mdful;
    for(j in 1:nln)  for(k in 1:(nfrg+1))  for(i in 1:(nmass+1-k)) tmp[j,i+k-1]<-tmp[j,i+k-1]-corr[i]*(fr[j,k]);
      fr<-mdistr(nfrg,tmp,mmteor,nln);
        }}
  if(md=="co") {# SUMMED INJECTIONS:
   for(j in 1:nlnsum)  gcmsum[j,1:nmass]<-(gcmsum[j,1:nmass]-corr)
     frsum<-mdistr(nfrg,gcmsum,mmteor,nlnsum)
     # ALL INJECTIONS:
                 for(j in 1:nln)  mdful[j,1:nmass]<-mdful[j,1:nmass]-corr;
     fr<-mdistr(nfrg,mdful,mmteor,nln)
               }
write("*** All samples fully corrected **",fn1,append=T)
  write.table(cbind(id,round(fr,5)),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F)
        write("\n",fn1,append=TRUE);
#
     i<-1 # Statistics for SUMMED INJECTIONS
  while(i<nlnsum){ cnlab<-substr(idsum[i],1,6); k<-i
    while(substr(idsum[k+1],1,6)==cnlab) k<-k+1
    # mean for SUMMED INJECTIONS
      if(k>i) {sredn<-frsum[1,]; sid<-paste("mean",substr(idsum[i],1,6),sep="_")
                 for(j in 1:nmass) sredn[j]<-mean(frsum[i:k,j])
                 for(j in 1:nmass) if(sredn[j]<0.) sredn[j]<-0.
                 srsum<-sum(sredn[1:nmass])
               sredn[1:nmass]<-sredn[1:nmass]/srsum # renormalization
    # SD for SUMMED INJECTIONS
               otkl<-frsum[1,]; oid<-paste("sd",substr(idsum[i],1,6),sep="_");
                 for(j in 1:nmass) otkl[j]<-sd(frsum[i:k,j]);
        frsum<-rbind(frsum[1:k,],sredn,otkl,frsum[(1+k):nrow(frsum),]);
        idsum<-c(idsum[1:k],sid,oid,idsum[(1+k):length(idsum)]);
    write.table(t(c(sid,nfrg,round(sredn[1:(nfrg+1)],5))),fn2,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
    write.table(t(c(oid,nfrg,round(otkl[1:(nfrg+1)],5))),fn2,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
# write.table(t(c("*Correction: ",round(corr,5))),fn2,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
                  i<-k+2;}
          i<-i+1; nlnsum<-nrow(frsum);}
# write data:
 write("*** Statistics, samples fully corrected **",fn1,append=TRUE)
    write.table(cbind(idsum,round(frsum,5)),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
 write.table(t(c("*Correction: ",round(corr,5))),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
        write("\n",fn1,append=TRUE);
         return(list(fr,nfrg,iln))
}

findfrg<-function(rada,iln,colmet,colfrg){
 i<-iln; k<-i; l<-1; chast<-rada[1,]; len<-nrow(rada);
  metab<-rada[i,colmet]; fra<-rada[i,colfrg]; fra1<-fra;
  while(metab==rada[i,colmet]){ #count of lines for: i-> metabolite, k->fragments
     if(fra==rada[i,colfrg]) {k=k+1; chast=rbind(chast,rada[i,])} #bind lines for fragment
       else{ iii=0; for(ii in 1:l) if(rada[i,colfrg] != fra1[ii]) iii=iii+1
              if(iii==l){l=l+1; fra1[l]=rada[i,colfrg]}
              }
                           if(i==len) break;  i=i+1; }
   return(list(chast,fra1,i,k)) }

exfrag<-function(rada,frag,iln,colmet,colfrg){
   i=iln; chast<-rada[1,]; len=length(rada[,1]);
   while(i<=len) {
     if(frag==rada[i,colfrg]){ chast=rbind(chast,rada[i,]); i=i+1 }
      else { i=i+1}
     }
   return(chast)
}
run_midcor<-function(inputFileName="xramidout.csv", output="midcorout.csv",mode="con"){
  fn<-file.path(inputFileName);
  fn1<-paste(fn,"_c",sep="");
  fn2<-file.path("edata");
  write("",fn1);
  write("",output);
  rada<-read.table(fn, sep=",");   # read experimental data
  for(i in 1:ncol(rada)) {
        if(grepl("time",rada[1,i])) {coltime=i}  # column of incubation time
         if(grepl("Metab",rada[1,i])) {colmet=i}  # column of metabolite name
        if(grepl("atomic pos",rada[1,i])) {colfrg=i}  # carbon positions in the fragment
        if(grepl("signal intensity",rada[1,i])) {inten=i}  #  column of signal intensity
        if(grepl("isotopol",rada[1,i])) {isoname=i; } # isotopolog symbols
        if(grepl("abundance",rada[1,i])) {abund=i}  # calculated fractions of isotopologs
  }
  tt<-na.omit(as.numeric(levels(rada[,coltime])))[1]
  cat(c(tt,-1),"\n",c(1,2),"growth_0_t\nAnusha\n",file=fn2);
   iln=2;
tot=data.frame();
   Mtb<-levels(rada[,colmet]); numet<-length(Mtb)-1
   for(ii in 1:numet){
        write((paste("----",Mtb[ii]," ---")),fn1,append=TRUE)
        a=findfrg(rada,iln,colmet,colfrg); # a: 1-part of rada referred to selected fragment
# 2-fragments, 3-number of lines for the metabolite, 4-for fragment
        chast=a[[1]]; fra1=a[[2]]; i=a[[3]]; k=a[[4]]; lnfrg<-length(fra1)
        if(lnfrg>1) lnst=iln
      for(frcyc in 1:lnfrg){
        if(frcyc>1) chast=exfrag(rada,fra1[frcyc],lnst,colmet,colfrg);
        res<-correct(chast,fn1,fn2,mode)#,"var")
        nstrok<-nrow(res[[1]])
        i=1; newcol<-numeric()
              for(j in 1:nstrok){
                while(!grepl("C0",chast[i,isoname])) {i=i+1;}
                        for(k in 1:(res[[2]]+1)) {
                                  newcol[i]<-res[[1]][j,k]*100.; i=i+1;
                                                     }
              }
           newcol<-round(newcol,3);   if(length(newcol)<nrow(chast)) newcol<-c(newcol,0)
        chast=cbind(chast[-c(1),],newcol[-c(1)]); iln<-iln+res[[3]]-2; tot=rbind(tot,chast)
        }
    }
    alname=rada[1,];
#    alname=cbind(rada[1,],"midcor")
#    tot[,abund+1]=as.numeric(as.character(tot[,abund+1]))
    tot<-tot[-(abund)]
#    tot[,inten]=as.numeric(as.character(tot[,inten]))
    tot[,inten]=as.numeric(levels(tot[,inten])[tot[,inten]])
   write.table(alname,output,sep=",",append=F,col.names=FALSE, row.names = F);
   write.table(tot,output,sep=",",append=TRUE,col.names=FALSE, row.names = F);
   return(newcol) }


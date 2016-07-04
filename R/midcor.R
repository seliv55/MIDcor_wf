correct<-function(chast,fn1,fn2,mdcor="con"){
  md<-substr(mdcor,1,2);
  onemet=convert(chast,2); #datacont=c("id","mm","i","nmet","nC","nfrg")
  id=onemet[[1]]; mm=onemet[[2]]; iln=onemet[[3]]; nmet<-onemet[[4]]; nC<-onemet[[5]];
  nfrg<-onemet[[6]]; mmsum=onemet[[7]]; idsum=onemet[[8]]; nSi<-0; nS<-0;
   ncol=length(mm[1,]);  nmass=ncol-1; nln<-length(id);
# normalization
   mdf=data.frame(cbind(id,mm), row.names = NULL);
#    mdf[,1]=as.character(mdf[,1]);
    for(i in 2:(1+ncol)) mdf[,i]=as.numeric(as.character(mdf[,i]));
    ef<-sum(mdf[2])/sum(mdf[3]);
         gcmsn<-norm(mdf,ef); mdful=gcmsn; nstrok=length(mdful[,1]);
# theoretic distribution:
          mmteor<-mtr(nfrg,nmass,nC,nSi,nS);
# mass fractions
   fr<-mdistr(nfrg,gcmsn,mmteor,nln);# write mass fractions without correction:
 write("*** MID for each injection, corrected only for natural 13C, 29,30Si, 33,34S ***",fn1,append=T)
  write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
        write("\n",fn1,append=TRUE);

  mdf=data.frame(cbind(idsum,mmsum), row.names = NULL);
  mdf[,1]=as.character(mdf[,1]);
  for(i in 2:(1+ncol)) mdf[,i]=as.numeric(as.character(mdf[,i]));
# normalization
         gcmsn<-norm(mdf,ef); gcmsn1<-gcmsn; nln=length(mdf[,1]);
# mass fractions
      fr<-mdistr(nfrg,gcmsn,mmteor);
 write("*** Summed injections for each plate, corrected only for natural 13C, 29,30Si, 33,34S **",fn1,append=T)
  write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
        write("\n",fn1,append=TRUE);
         gcmsn<-gcmsn1;
# correction
     corr<-numeric(nmass);
     for(i in 1:nln) {if(grepl("cold_c",mdf[i,1])) break;}
       corr<-gcmsn[i,3:(2+nmass)]-mmteor[1,];
 if(md=="va") { for(ii in 1:9) {tmp<-gcmsn;
    for(j in 1:nln)  for(k in 1:(nfrg+1))  for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k]);
      fr<-mdistr(nfrg,tmp,mmteor,nln);
    for(j in 1:nstrok)  for(k in 1:(nfrg+1))  for(i in 1:(nmass+1-k)) mdful[j,i+1+k]<-mdful[j,i+1+k]-corr[i]*(fr[j,1+k]);
      mdful<-mdistr(nfrg,mdful,mmteor,nstrok);
        }}
  if(md=="co") { for(j in 1:nln)  gcmsn[j,3:(2+nmass)]<-gcmsn[j,3:(2+nmass)]-corr;
     fr<-mdistr(nfrg,gcmsn,mmteor,nln);
                 for(j in 1:nstrok)  mdful[j,3:(2+nmass)]<-mdful[j,3:(2+nmass)]-corr;
     mdful<-mdistr(nfrg,mdful,mmteor,nstrok);}

                             i<-1; lfr<-length(fr);
  while(i<nln){ cnlab<-as.character(fr[i,1]); k<-i;
    while(as.character(fr[k+1,1])==cnlab) k<-k+1;
      if(k>i) {sredn<-fr[1,]; sredn[1]<-paste(chast[2,nmet],chast[2,nmet+2],fr[i,1],"mean",sep="_");
                 for(j in 2:lfr) sredn[j]<-mean(fr[i:k,j]);
                 for(j in 2:lfr) if(sredn[j]<0.) sredn[j]<-0.;
                 srsum<-sum(sredn[2:(2+nfrg)]);
               sredn[2:(2+nfrg)]<-sredn[2:(2+nfrg)]/srsum;
               otkl<-fr[1,]; otkl[1]<-paste(chast[2,nmet],chast[2,nmet+2],fr[i,1],"sd",sep="_");
                 for(j in 2:length(otkl)) otkl[j]<-sd(fr[i:k,j]);
        fr<-rbind(fr[1:k,],sredn,otkl,fr[(1+k):nrow(fr),]);
    write.table(format(sredn,digits=4),fn2,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
    write.table(format(otkl,digits=4),fn2,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
    cor1=cbind("*Correction: ",corr)
 write.table(format(cor1,digits=4),fn2,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
                  i<-k+2;}
          i<-i+1; nln<-nrow(fr);}
# write data:
       razn<-gcmsn[1,];  razn[1]<-"correction"; for(i in 1:nmass) razn[i+1]<-corr[i];
       frr<-fr[1:(length(fr)-2)]
 write("*** Statistics, samples fully corrected **",fn1,append=TRUE)
    write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
 write("*** Correction factor: **",fn1,append=TRUE)
     write.table(format(razn,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
        write("\n",fn1,append=TRUE);
         return(list(mdful,nfrg,iln))
}

findfrg<-function(rada,iln,colmet,colfrg){
 i=iln; k=i; l=1; chast<-rada[1,]; len=length(rada[,1]);
  metab=rada[i,colmet]; fra=rada[i,colfrg]; fra1=fra;
  while(metab==rada[i,colmet]){
     if(fra==rada[i,colfrg]) {k=k+1; chast=rbind(chast,rada[i,])}
       else{ iii=0; for(ii in 1:l) {if(rada[i,colfrg] != fra1[ii]) {iii=iii+1}}
              if(iii==l){l=l+1; fra1[l]=rada[i,colfrg]}}
                           if(i==len){break;};  i=i+1; }
   return(list(chast,fra1,i,k)) }

exfrag<-function(rada,frag,iln,colmet,colfrg){
   i=iln; chast<-rada[1,]; len=length(rada[,1]);
   while(i<=len) {
     if(frag==rada[i,colfrg]){ chast=rbind(chast,rada[i,]); i=i+1 }
      else { i=i+1}
     }
   return(chast)
}
run_midcor<-function(inputFileName, output){
  fn<-file.path(inputFileName);
  fn1<-paste(fn,"_c",sep="");
  fn2<-file.path("isodyn");
  write("",fn2);
  write("",fn1);
  write("",output);
  rada<-read.table(fn, sep=",");   # read experimental data
  for(i in 1:length(rada)) {
        if(grepl("isotop",rada[1,i])) {isoname=i; } # column of signal intensity
        if(grepl("Metab",rada[1,i])) {colmet=i}
        if(grepl("atomic pos",rada[1,i])) {colfrg=i}
  }
   iln=2;
tot=data.frame();
   Mtb=levels(rada[,colmet]); numet=length(Mtb)-1;
   for(ii in 1:numet){
        write(as.character(paste("----",Mtb[ii]," ---")),fn1,append=TRUE);
        a=findfrg(rada,iln,colmet,colfrg);
        chast=a[[1]]; fra1=a[[2]]; i=a[[3]]; k=a[[4]]; lnfrg=length(fra1);
        if(lnfrg>1) lnst=iln;
      for(frcyc in 1:lnfrg){
        if(frcyc>1) chast=exfrag(rada,fra1[frcyc],lnst,colmet,colfrg);
        res=correct(chast,fn1,fn2)#,"var")
        nstrok=length(res[[1]][,1]);
        i=1; newcol=as.character(chast[,isoname]);
              for(j in 1:nstrok){
                while(chast[i,isoname]!="m0") {i=i+1;}
                        for(k in 1:(res[[2]]+1)) {
                                  newcol[i]=res[[1]][j,k+1]; i=i+1;
                                                     }
              }
        chast=cbind(chast[-c(1),],newcol[-c(1)]); iln=iln+res[[3]]-2; tot=rbind(tot,chast)
        }
    }
    alname=cbind(rada[1,],"midcor")
   write.table(alname,output,sep=",",append=F,col.names=FALSE, row.names = F);
   write.table(tot,output,sep=",",append=TRUE,col.names=FALSE, row.names = F);
   return(newcol) }


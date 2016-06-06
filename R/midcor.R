correct<-function(data,fn1,mdcor="con"){
  md<-substr(mdcor,1,2);
  id=data[[1]]; mm=data[[2]]; iln=data[[3]];
   nmet<-data[[4]]; nC<-data[[5]]; nfrg<-data[[6]]; nSi<-0; nS<-0;
   ncol=length(mm[1,]);  nmass=ncol-1; nln<-length(id);
# normalization
    mdf=data.frame(cbind(id,mm), row.names = NULL);
    mdf[,1]=as.character(mdf[,1]);
    for(i in 2:(1+ncol)) mdf[,i]=as.numeric(as.character(mdf[,i]));
    ef<-sum(mdf[2])/sum(mdf[3]);
         gcmsn<-norm(mdf,ef); mdful=gcmsn; nstrok=length(mdful[,1]);
# theoretic distribution:
          mmteor<-mtr(nfrg,nmass,nC,nSi,nS);
# mass fractions
   fr<-mdistr(nfrg,gcmsn,mmteor,nln);# write mass fractions without correction:
 write("\n\n*** MID for each injection, corrected only for natural 13C, 29,30Si, 33,34S ***",fn1,append=TRUE)
  write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
# Sum the data for various injections from the same plate:
               i<-1; nln=length(mdf[,1]);
  while(i<(nln)){ l=nchar(mdf[i,1]); cnlab<-substr(mdf[i,1],l-5,l-3);
      while(substr(mdf[i+1,1],l-5,l-3)==cnlab) {
      mdf[i,2:(ncol+1)]<-mdf[i,2:(ncol+1)]+mdf[i+1,2:(ncol+1)];
        mdf<-mdf[-(i+1),]; nln<-nln-1; if(i==nln){break;}}
             i<-i+1;}

# normalization
         gcmsn<-norm(mdf,ef); gcmsn1<-gcmsn; nln=length(mdf[,1]);
# mass fractions
      fr<-mdistr(nfrg,gcmsn,mmteor);
 write("\n\n*** Summed injections for each plate, corrected only for natural 13C, 29,30Si, 33,34S **",fn1,append=TRUE)
  write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
         gcmsn<-gcmsn1;
          mmteor<-mtr(nfrg,nmass,nC,nSi,nS);
# correction
     corr<-numeric(nmass);
     for(i in 1:nln) {if(grepl("coldc",mdf[i,1])) break;}
       corr<-gcmsn[i,3:(2+nmass)]-mmteor[1,];
 if(md=="va") { for(ii in 1:9) {tmp<-gcmsn;
    for(j in 1:nln)  for(k in 1:(nfrg+1))  for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k]);
      fr<-mdistr(nfrg,tmp,mmteor,nln);
        }}
  if(md=="co") { for(j in 1:nln)  gcmsn[j,3:(2+nmass)]<-gcmsn[j,3:(2+nmass)]-corr;
     fr<-mdistr(nfrg,gcmsn,mmteor,nln);
      for(j in 1:nstrok)  mdful[j,3:(2+nmass)]<-mdful[j,3:(2+nmass)]-corr;
     mdful<-mdistr(nfrg,mdful,mmteor,nstrok);}

#                             i<-1; lfr<-length(fr);
#  while(i<nln){ cnlab<-substr(as.character(fr[i,1]),4,5); k<-i;
#    while(substr(as.character(fr[k+1,1]),4,5)==cnlab) k<-k+1;
#      if(k>i) {sredn<-fr[1,]; sredn[1]<-as.factor("**mean**");
#                 for(j in 2:lfr) sredn[j]<-mean(fr[i:k,j]);
#                 for(j in 2:lfr) if(sredn[j]<0.) sredn[j]<-0.;
#                 srsum<-sum(sredn[2:(2+nfrg)]);
#               sredn[2:(2+nfrg)]<-sredn[2:(2+nfrg)]/srsum;
#               otkl<-fr[1,]; otkl[1]<-as.factor("**sd**");
#                 for(j in 2:length(otkl)) otkl[j]<-sd(fr[i:k,j]);
#        fr<-rbind(fr[1:k,],sredn,otkl,fr[(1+k):nrow(fr),]);
#                  i<-k+2;}
#          i<-i+1; nln<-nrow(fr);}
# write data:
       razn<-gcmsn[1,];  razn[1]<-"correction"; for(i in 1:nmass) razn[i+1]<-corr[i];
       frr<-fr[1:(length(fr)-2)]
 write("\n\n*** Statistics, samples fully corrected **",fn1,append=TRUE)
    write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
 write("*** Correction factor: **",fn1,append=TRUE)
     write.table(format(razn,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
         return(mdful)
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
# read experimental data
run_midcor<-function(inputFileName, output){
  fn<-file.path(inputFileName);
  fn1<-paste(fn,"_c",sep="");
  write("",fn1);
  rada<-read.table(fn, sep=",");
  for(i in 1:length(rada)) {
        if(grepl("isotop",rada[1,i])) {isoname=i; } # column of signal intensity
        if(grepl("Metab",rada[1,i])) {colmet=i}
        if(grepl("atomic pos",rada[1,i])) {colfrg=i}
  }
   iln=2;
tot=data.frame();
   Mtb=levels(rada[,colmet]); numet=length(Mtb)-1;
   for(ii in 1:numet){
        write(as.character(paste("\n----",Mtb[ii]," ---\n")),fn1,append=TRUE);
        a=findfrg(rada,iln,colmet,colfrg);
        chast=a[[1]]; fra1=a[[2]]; i=a[[3]]; k=a[[4]]; lnfrg=length(fra1);
        if(lnfrg>1) lnst=iln;
      for(frcyc in 1:lnfrg){
        if(frcyc>1) chast=exfrag(rada,fra1[frcyc],lnst,colmet,colfrg);

        data=convert(chast,2); #datacont=c("id","mm","i","nmet","nC","nfrg")
        res=correct(data,fn1)
        nstrok=length(res[,1]);
        i=1; newcol=as.character(chast[,isoname]);
        for(j in 1:nstrok){
                while(chast[i,isoname]!="m0") {i=i+1;}
                        for(k in 1:(data[[6]]+1)) {
                                  newcol[i]=as.character(res[j,k+1]); i=i+1;
                        }
        }
        chast=cbind(chast,newcol); iln=iln+data[[3]]-2; tot=rbind(tot,chast)
        }
    }
     write.table(tot,output,sep=",",append=TRUE,col.names=FALSE, row.names = F);
 return(newcol) }


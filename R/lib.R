binom<-function(n){ #calculation of binomial cefficients
 bf<-numeric(n); bf[1]<-1;
 if(n>1){ for(i in 2:n) bf[i]<-bf[i-1]*i;  bc<-numeric(n);
   for(i in 1:(n-1)) bc[i]=bf[n]/bf[i]/bf[n-i];}
   else  bc<-1;
     bc }

carbons<-function(n, bc) {#natural distribution of 13C according to the number of carbons in whole fragment
  pc12<-0.989; pc13<-0.011;
    mc<-numeric(n+1); mc[1]<-pc12^n;
 if(n>1) {for(i in 1:(n-1)) mc[i+1]<-bc[i]*pc12^(n-i)*pc13^i; mc[n+1]<-pc13^n;}
   else mc<-1;
       mc }

silicium<-function(mc){#correction of natural mass distribution accounting for Si
  pSi0<-0.9223; pSi1<-0.0467; pSi2<-0.031;
   n<-length(mc)
    m<-numeric(n+2)
     m[1]<-mc[1]*pSi0;
     m[2]<-mc[2]*pSi0+mc[1]*pSi1;
  for(i in 3:n) m[i]<-mc[i]*pSi0+mc[i-1]*pSi1+mc[i-2]*pSi2;
     m[n+1]<-mc[n]*pSi1+mc[n-1]*pSi2;
       m[n+2]<-mc[n]*pSi2;
         m }

sulfur<-function(mc){#correction of natural mass distribution accounting for S
  pS0<-0.9504; pS1<-0.0075; pS2<-0.0421;
   n<-length(mc)
    m<-numeric(n+2)
     m[1]<-mc[1]*pS0;
     m[2]<-mc[2]*pS0+mc[1]*pS1;
  for(i in 3:n) m[i]<-mc[i]*pS0+mc[i-1]*pS1+mc[i-2]*pS2;
     m[n+1]<-mc[n]*pS1+mc[n-1]*pS2;
       m[n+2]<-mc[n]*pS2;
         m }

tdist<-function(nc,nsi,ns){ #final natural mass distribution
   bcc<-binom(nc);
    m<-carbons(nc, bcc);
 if(nsi>0) for(i in 1:nsi)  m<-silicium(m);
 if(ns>0) for(i in 1:ns)  m<-sulfur(m);
      m }

norm<-function(mmsd,f1){#normalization of GCMS data accounting for the loss of protons
        msd<-mmsd; ncol<-length(msd);# f1<-f1[1];
    for (i in 3:(ncol-1)) msd[i]<-(msd[i]+msd[i]*f1-msd[i+1]*f1);
   scol<- msd[3]; for(i in 4:(ncol-1)) scol<-scol+msd[i];
  msd[2]<-msd[2]*0;
 for (i in 3:(ncol-1)) msd[i]<-(msd[i]/scol); #normalization
  scol<- msd[3]; for(i in 4:(ncol-1)) scol<-scol+msd[i];#=1
    msd[ncol]<-scol;
      msd }

mtr<-function(nCfrg,nmass,nCder,nSi,nS){ #various numbers of labeled 13C with tails of natural distributions
    mt<-tdist(nCder,nSi,nS);
     lem<-length(mt); if(lem<nmass) for(i in (lem+1):nmass) mt[i]<-0.0;
      mt<-mt[1:nmass];
       mt<-mt/sum(mt);
        mm<-numeric(nmass*(nCfrg+1));
         attr(mm,"dim")<-c(nCfrg+1,nmass);
          mm[1,]<-mt;
       for(i in 1:nCfrg) { mt<-tdist(nCder-i,nSi,nS); lem<-length(mt);
         if(lem<(nmass-i)) {for(k in (lem+1):(nmass-i)) mt[k]<-0.0;}
             mt<-mt[1:(nmass-i)];  mt<-mt/sum(mt);
               for(j in 1:length(mt)) mm[i+1,j+i]<-mt[j];}
              mm}

mdistr<-function(nreal,msd,mm,nln){ #label incorporation
          fr<-msd; ncol<-length(fr); for(i in 2:ncol) fr[i]<-0;
 for(j in 1:(nreal+1)) {fr[j+1]<-msd[j+2]/mm[j,j];
   for(i in j:(ncol-3)) msd[i+2]<-msd[i+2]-fr[j+1]*mm[j,i];}
# normalization:
  fr[ncol]=fr[2];
        for(i in 3:(ncol-2)) {fr[ncol]<- fr[ncol]+fr[i];}

     for(i in 2:(ncol-2)) fr[i]<- fr[i]/fr[ncol];
      fr }

ginfo<-function(rawdat,ln){
  for(i in 1:length(rawdat)){if(grepl("deriv",rawdat[1,i])) {kk=i;} # column of formula
              else {if(grepl("posit",rawdat[1,i])) {ifrg=i;}} #  column of studied fragment
                           }
       a=strsplit(as.character(rawdat[ln,kk]),"C")[[1]][2]; # C atoms in derivate
       nCder=as.numeric(strsplit(a,"H")[[1]][1]);

       a=strsplit(as.character(rawdat[ln,ifrg]),"C")[[1]]; # C atoms in the fragment
       nCfrg=as.numeric(a[3])-as.numeric(strsplit(a[2],"-")[[1]][1])+1;
       return (list(nCder,nCfrg,ifrg))
       }
convert<-function(rdat,iln){
        colid=1; a=ginfo(rdat,iln);
         nCder=a[[1]]; nCfrg=a[[2]]; ifrg=a[[3]];
  for(i in 1:length(rdat)){if(grepl("intens",rdat[1,i])) {mdis=i;} # column of signal intensity
               else {if(grepl("Metab",rdat[1,i])) nmet=i;} #  column of studied fragment
         if (grepl("Inject",rdat[1,i])) {inj=i; repl=i+1; cond=i+2;}
         }
         i=iln;
         id=as.character(rdat[i,colid]);
            samp=rdat[i,inj];
             met=rdat[i,nmet];
             frag=rdat[i,ifrg];
          m=numeric(); k=1;
       while(samp==rdat[i,inj]) {if(frag==rdat[i,ifrg]){ m[k]=as.numeric(as.character(rdat[i,mdis])); i=i+1; k=k+1;}
             else {i=i+1;}};
       lm=length(m);
       mm=matrix(nrow=1,ncol=lm); mm[1,]=m;
       while(met==rdat[i,nmet]){ id=c(id,as.character(rdat[i,colid]));
         samp=rdat[i,inj];
            m=numeric(); k=1;
       while(samp==rdat[i,inj]) {
         {
           if(frag==rdat[i,ifrg]){
             m[k]=as.numeric(as.character(rdat[i,mdis]));
             i=i+1;
             k=k+1;
           }
           else {
             i=i+1
           }
         }
       if(i>length(rdat[,1])) {break}}
       if(lm==length(m)) {mm=rbind(mm,m)};
       if(i>length(rdat[,1])) {break}}
       return(list(id,mm,i,nmet,nCder,nCfrg))
}



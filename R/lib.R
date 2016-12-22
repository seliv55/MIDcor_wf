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

Norm<-function(msd,f1){#normalization of GCMS data accounting for the loss of protons
        nucol<-ncol(msd);# f1<-f1[1];
    for (i in 2:(nucol-1)) msd[,i]<-(msd[,i]+msd[,i]*f1-msd[,i+1]*f1);
    msd<-msd[,-1]; nucol<-nucol-1
   scol<- msd[,1]; for(i in 2:(nucol-1)) scol<-scol+msd[,i];
 for (i in 1:(nucol-1)) msd[,i]<-(msd[,i]/scol); #normalization
  scol<- msd[,1]; for(i in 2:(nucol-1)) scol<-scol+msd[,i];#=1
    msd[,nucol]<-scol;
      return(msd) }

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
          fr<-msd; nucol<-ncol(fr); for(i in 1:nucol) fr[,i]<-0;
 for(j in 1:(nreal+1)) {fr[,j]<-msd[,j]/mm[j,j];
   for(i in j:(nucol-1)) msd[,i]<-msd[,i]-fr[,j]*mm[j,i];}
# normalization:
  fr[,nucol]=fr[,1];
        for(i in 2:(nucol-1)) {fr[,nucol]<- fr[,nucol]+fr[,i];}

     for(i in 1:(nucol-1)) fr[,i]<- fr[,i]/fr[,nucol];
  fr[,nucol]=fr[,1];
        for(i in 2:(nucol-1)) {fr[,nucol]<- fr[,nucol]+fr[,i];}
      return(fr) }

ginfo<-function(rawdat,ln){
  for(i in 1:ncol(rawdat)){if(grepl("formula derivatized",rawdat[1,i])) {kk=i;} # column of formula
              else {if(grepl("atomic positions",rawdat[1,i])) {colfrg=i;}} #  column of studied fragment
                           }
       a=strsplit(as.character(rawdat[ln,kk]),"C")[[1]][2]; # C atoms in derivate
       nCder=as.numeric(strsplit(a,"H")[[1]][1]);

       a=strsplit(as.character(rawdat[ln,colfrg]),"C")[[1]]; # C atoms in the fragment
       nCfrg=as.numeric(a[3])-as.numeric(strsplit(a[2],"-")[[1]][1])+1;
       return (list(nCder,nCfrg,colfrg))  }
       
convert<-function(rdat,iln){
        colid=1; a=ginfo(rdat,iln);colmet=0;
         nCder=a[[1]]; nCfrg=a[[2]]; colfrg=a[[3]];
  for(i in 1:ncol(rdat)){if(grepl("signal intens",rdat[1,i])) {coldis=i} # column of signal intensity
               else {if(grepl("Metab",rdat[1,i])) colmet=i;} #  column of metabolite name
         if (grepl("Inject",rdat[1,i])) {colinj=i; colrep=i+1;}
         if (grepl("labelled pos",rdat[1,i])) colab=i
         }
         ninj<-rdat[iln,colinj];
         nrepl<-rdat[iln,colrep];
         met<-rdat[iln,colmet];
         frag<-rdat[iln,colfrg];
         m<-numeric();     k<-1;
           i<-iln;  id<-paste(rdat[i,colid],substr(rdat[i,colab],1,3),sep="_");   idsum<-id
       while(ninj==rdat[i,colinj]) {# distribution of intensities for 1st injection
         if(frag==rdat[i,colfrg]){
               m[k]=as.numeric(as.character(rdat[i,coldis])); i=i+1; k=k+1}
         else i=i+1 }
       lm<-length(m)
       mm<-matrix(nrow=1,ncol=lm); mmm=mm # mm - for all inj; mmm - for summed inj
        mm[1,]=m; m1=m                    # m - for all inj; m1 - for summed inj
       while(met==rdat[i,colmet]){  #  metabolite cycle
        id<-c(id,as.character(rdat[i,colid]))
         ninj<-rdat[i,colinj];
            m<-numeric(lm); k<-1
             if(nrepl!=rdat[i,colrep]){ newrep<-paste(rdat[i,colid],substr(rdat[i,colab],1,3),sep="_")
              idsum=c(idsum,newrep);
               mmm=rbind(mmm,m1); nrepl=rdat[i,colrep]; m1=m;}
       while(ninj==rdat[i,colinj]) {  #  replicate cycle
          if(frag==rdat[i,colfrg]){
             m[k]=as.numeric(as.character(rdat[i,coldis]));
              if(nrepl==rdat[i,colrep]) m1[k]=m1[k]+m[k] # sum of injections
             i=i+1;  k=k+1 }
           else { i=i+1 }
       if(i>nrow(rdat)) {mmm=rbind(mmm,m1); break}}
       if(lm==length(m)) {mm=rbind(mm,m)};
       if(i>nrow(rdat)) {break}}
       return(list(id,mm,i,colmet,nCder,nCfrg,mmm[-1,],idsum))
}

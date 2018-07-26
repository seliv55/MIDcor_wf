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

ginfo<-function(rawdat,ifor,colfrg){
       a<-as.character(rawdat[2,ifor]) 
       pC<-regexpr("C",a)  # C-position in formula
       pH<-regexpr("H",a)  # H-position in formula
       nCder<-as.numeric(substr(a,pC+1,pH-1)); # C atoms in derivate
       nSi<-0; nS<-0
       pSi<-regexpr("Si",a); if(pSi>0) nSi<-as.numeric(substr(a,pSi+2,pSi+2)) # Si atoms in derivate
       pS<-gregexpr("S",a)[[1]]; if(length(pS)>1) { # S atoms in derivate
          if(pS[1]!=pSi) nS<-as.numeric(substr(a,pS[1]+1,pS[1]+1)) else nS<-as.numeric(substr(a,pS[2]+1,pS[2]+1))}
             
       a<-as.character(rawdat[2,colfrg])
       pCf<-gregexpr("C",a)[[1]]; # C-position in formula
       iCb<- as.numeric(substr(a,pCf[1]+1,pCf[1]+1))
       iCe<- as.numeric(substr(a,pCf[2]+1,pCf[2]+1))
       nCfrg<- iCe-iCb+1
return (list(nCder,nCfrg,nSi,nS))  }

convert<-function(rdat){
        colid<- 1;
  tit<-data.frame(lapply(rdat[1,], as.character), stringsAsFactors=FALSE)
        coldis<- grep("signal", tit)  # column of signal intensity
        colmet<- grep("Metab", tit)  # column of metabolite name
        coliso<- grep("isotopolog", tit)  # column of mass iso
        ifor<- grep("formula", tit)  # column of derivates
        colfrg<- grep("atomic", tit)  # column of fragments
        colab<- grep("labelled", tit)  # column of tracer
   a<- ginfo(rdat,ifor,colfrg);
         nCder<- a[[1]]; nCfrg<- a[[2]]; nSi<- a[[3]]; nS<- a[[4]]; 
  
         met<-rdat[2,colmet]
         frag<-rdat[2,colfrg]
        
           rdat[,coliso]<-as.character(rdat[,coliso])
           i<-2; id<-character(); first<-0;
            m00<-0; m<-numeric();  labmet<-character();
            
      while(i<nrow(rdat)){
      if((i<nrow(rdat))&(grepl("-1",rdat[i,coliso]))) {
      m00=as.numeric(as.character(rdat[i,coldis]))
       i<-i+1;    k<-1 }
       while(as.integer(strsplit(rdat[i,coliso],"13C")[[1]][2]) >= 0) {
        m[k]=as.numeric(as.character(rdat[i,coldis]));
         i<-i+1; if(!grepl("13C",rdat[i,coliso])) break; k<-k+1 }
          f<- m00/m[1]
          if(f<0.2) { for(eim in 1:(length(m)-1)) m[eim]<- (1+f)*m[eim]+f*m[eim+1]
             m<- m/sum(m);
             first<-first+1;
    if(first==1) { mm<-matrix(nrow=1,ncol=length(m)); mm[1,]<-m }
    else mm<-rbind(mm,m[1:ncol(mm)]);
      id<-c(id,as.character(rdat[i-3,colid]));
       labmet<-c(labmet,as.character(rdat[i-3,colab])) }
           else print(paste(rdat[3,colmet],' m00=',m00,' m1=',m[1]))
      while((i<nrow(rdat))&(!grepl("13C",rdat[i,coliso]))) i<-i+1
      }
           return(list(id,mm,nCder,nCfrg,nSi,nS,labmet))
}

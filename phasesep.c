/*
 model for predicting phase separation between negatively charged (RNA) particles 
 and positively charged (protein) particles

 Michael Feig, Bercem Dutagaci
 Michigan State University

 initial version April 2020

 reference:
 B. Dutagaci, G. Nawrocki, J. Goodluck, L. Lapidus, M. Feig:
 Charge-Driven Phase Separation of RNA and Proteins without Disorder
 bioRxiv(2020)
 

 compile with: 

 g++ -O3 -o phasesep phasesep.c

 usage: 

 phasesep [options]
    -verbose	 			// provide extra output
    -highres				// output with more digits
    -tag <name>				// add tag to output
    -temp <value>			// temperature in K
    -n <rnavalue> <posvalue> 		// number of molecules (can be fractional)
    -nrna <value> -npos <value>        
    -c <rnavalue> <posvalue> 		// concentrations in mM   
    -crna <value> -cpos <value>	
    -len <value>			// system size: len*len*len
    -kappa <value>			// salt screening; typical: 0.5-2.0 
     
    -q <rnavalue> <posvalue> 		// charges
    -qrna <value> -qpos <value>
    -r <rnavalue> <posvalue> 		// radii
    -rrna <value> -rpos <value>	
    -epsilon <value>			// model parameter; default: 4.0
    -a0 <value>				// model parameter; default: 3.0

    -rnarna <name>			// name for alternative RNA-RNA RDF profile
    -rnapos <name>			// name for alternative RNA-POS RDF profile
    -posrna <name>			// name for alternative POS-RNA RDF profile
    -pospos <name>			// name for alternative POS-POS RDF profile
    -rdfcut <value>			// cutoff for RDF beyond which it is set to 1
    -maxrad <value>			// radial integration limit 
    
    -thresh <value>			// numerical error tolerance; default: 0.0005
    -vfac <value>			// multiplier when scanning for volume; default: 0.9995
    -nrscale <value>			// step size when scanning for concentration; default: 0.02
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

int verbose=0;
int highres=0;

const double pi=3.14159265359;
const double kb=0.001987041;
const double navo=6.02214E23;
const double l10=log(10.0);

class Parameter {
  private:
    double sig(double r) { return 2.0*r/pow(2.0,1.0/6.0); } 
    double effcharge(double q) { return log((fabs(q)/fcharge)+1)*fcharge; }
    double spherevolume(double r) { return 4.0*pi/3.0*r*r*r; }    
    double asign(double q) { if (q>0.001 || q<-0.001) { return q/fabs(q); } else { return 1.0; } }
    double avalue(double q) { return asign(q)*sqrt((aa*fabs(effcharge(q)))); }
  
    double syslen;

    double rpos;
    double rrna;

  public:
    double temp;
    double kt2;

    double qpos;
    double qrna;
  
    double nrna;
    double npos;

    double kappa;
 
    double sysvol;

    double volpos;
    double volrna;

    double vfac;
    double nrscale;

    double fcharge;
    double aa;
    double a0;
    double epsilon;

    Parameter() : qpos(20.0), qrna(-75), nrna(297), npos(211), kappa(1.5), vfac(0.9999), nrscale(0.02), fcharge(20.0), aa(0.75), a0(3.0), epsilon(4.0) {
                  setsyslen(100.0); setrpos(3.5); setrrna(1.74); settemp(298.0); }

    void setsyslen(double slen) {
       syslen=slen;
       sysvol=slen*slen*slen;
    }

    void setrpos(double r) {
       rpos=r;
       volpos=spherevolume(r);
    }
 
    void setrrna(double r) {
       rrna=r;
       volrna=spherevolume(r);
    } 
    
    void settemp(double t) {
       temp=t;
       kt2=t*kb*2.0;
    }
 
    double rnadens() { return (double)nrna/sysvol; }
    double posdens() { return (double)npos/sysvol; }

    double qeffrna() { return effcharge(qrna); }
    double qeffpos() { return effcharge(qpos); }
    
    double sigrna() { return sig(rrna); }
    double sigpos() { return sig(rpos); }
 
    double dsigrna() { return sig(1.74); }
    double dsigpos() { return sig(3.5); }
 
    double arna() { return avalue(qrna); }
    double apos() { return avalue(qpos); }
};

class Energy {
  private:
    double sigma;
    double a;
    Parameter p;

  public:
    Energy() : sigma(1.0), a(1.0) {};
    Energy(double sigval, double aval, Parameter &par) : sigma(sigval), a(aval), p(par) {};

    double value(double r) { 
      return (4.0*p.epsilon*(pow(sigma/r,10)-pow(sigma/r,5))+(a+p.a0)*p.kappa*sigma/r*exp(-r/(p.kappa*sigma)))*0.239; 
    }
};    

class Field {
 private:
  static const int maxfields;

  char *tstr;
  char *fstr[300];

  int nfields;

  int first;

  int contains(char c, char *set);

 public:
  Field(char *str, char *sep=(char *)" \t"); 
  ~Field();

  void shift();


  int number() { return nfields-first; }

  char *operator[](int inx);
};

// ##### Field ######################################################

const int Field::maxfields=200;

int Field::contains(char c, char *set) {
  for (char *iptr=set; *iptr!=0; iptr++) {
    if (c==*iptr) {
      return 1;
    }
  }
  return 0;
}


Field::Field(char *inp, char *sep) : nfields(0), first(0) {
  tstr=new char[strlen(inp)+2];

  strcpy(tstr,inp);

  char *cptr=tstr;
  do {
    for (; *cptr!=0 && contains(*cptr,sep); cptr++);
    if (*cptr!=0) {
      fstr[nfields++]=cptr;
      for (;*cptr!=0 && !contains(*cptr,sep); cptr++);
      if (*cptr!=0) {
        *cptr=0;
        cptr++;
      }
    }
  } while (*cptr!=0);
}


Field::~Field() {
  delete tstr;
}

void Field::shift() {
  first++;
}

char *Field::operator[](int inx) {
  return fstr[inx-first];
}

// ##### RDF #############################################################

enum RDFType { DISPERSE=0, DILUTE=1, CLUSTER=2 };

class RDFData {
  private:
    int max;
    int n;
    double *r;
    double *g;

    double integralValue;

    double lookup(double rad, double rsig); 

  public:
    RDFData(int maxn=1000) : max(maxn), n(0) { 
      r=new double[maxn];
      g=new double[maxn];
    }
    ~RDFData() {
      delete r;
      delete g;
    } 
    
    void read(char *fname, double cut);
    void integrate(double sig, double refsig, double a, Parameter &par, double max);
    
    double precalcIntegral() { return integralValue; }
}; 

void RDFData::read(char *fname, double cut) {
  FILE *fptr;
  fptr=fopen(fname,"r");
  
  if (fptr==0) {
    fprintf(stderr,"cannot open RDF file %s\n",fname);
    exit(1);
  }

  double rval;
  double gval; 

  double lastval=-1.0;
  double cutg=-1.0;

  while (!feof(fptr)) {
    if (fscanf(fptr,"%lf%lf",&rval,&gval)) {
      if (n<max) {
        r[n]=rval/10.0;
        g[n]=gval;
        n++;

        if (rval>cut*10.0 && cutg<0) {
          cutg=lastval;
        }
        lastval=gval; 
      }      
    }
  }
  fclose(fptr);

  for(int i=0; i<n; i++) {
    if (r[i]>cut) {
      g[i]=1.0;
    } else {
      g[i]=g[i]/cutg;
    }
  }
}   

double RDFData::lookup(double rad, double rsig) {
  if (rad<=r[0]*rsig) { return 0.0; }
  if (rad>=r[n-1]*rsig) { return g[n-1]; }

  double val=-1.0;
  for (int i=0; i<n-1 && val<0; i++) {
    if (rad>=r[i]*rsig && rad<r[i+1]*rsig) {
      double drdf=g[i+1]-g[i];
      double dr=r[i+1]*rsig-r[i]*rsig;
      double xr=(rad-r[i]*rsig)/dr;
      val=xr*drdf+g[i];
    }     
  }

  if (val<0) {
    val=g[n-1];
  }

  return val;
}

void RDFData::integrate(double sig, double refsig, double a, Parameter &par, double max) {
  Energy ener(sig,a,par);

  int irfrom=2;
  int irto=int(max*10);
  
  int nint=int((irto-irfrom+2)/2);
  
  double *x=new double[nint];
  double *v=new double[nint];

  for (int i=0; i<nint; i++) {
    int ir=irfrom+i*2;
    double rad=(double)ir*0.1;
    x[i]=rad;
    v[i]=lookup(rad,sig/refsig)*ener.value(rad)*rad*rad;
  }
  
  double sum=0.0;
  double deltax=(x[nint-1]-x[0])/(double)(nint-1);
  
  sum+=v[0];
  for (int i=1; i<nint-1; i+=2) {
    sum+=4.0*v[i];
  }
  for (int i=2; i<nint-1; i+=2) {
    sum+=2.0*v[i];
  }
  sum+=v[nint-1];
   
  integralValue=4.0*pi*sum*deltax/3.0;   

  delete x;
  delete v;
}   

class RDF {
  private:
    RDFData *rdf;

  public:
    RDF(char *fname, double sig, double refsig, double a, Parameter &par, double cutval, double maxval) {
      rdf=new RDFData[3];

      char tname[512];

      sprintf(tname,"%s.c.rdf",fname);
      rdf[CLUSTER].read(tname,cutval);
      rdf[CLUSTER].integrate(sig,refsig,a,par,maxval);

      sprintf(tname,"%s.d.rdf",fname);
      rdf[DISPERSE].read(tname,cutval);
      rdf[DISPERSE].integrate(sig,refsig,a,par,maxval);

      sprintf(tname,"%s.o.rdf",fname);
      rdf[DILUTE].read(tname,cutval);
      rdf[DILUTE].integrate(sig,refsig,a,par,maxval);
    }       
    
    ~RDF() {
      delete[] rdf;
    }

    double pint(RDFType type) { return rdf[type].precalcIntegral(); }
};

class EnergyDens {
  private:
    double svol(double vratio, double temp);

  public:
    double xrr;
    double xrp; 
    double xpr;
    double xpp;
  
    EnergyDens() : xrr(0.0), xrp(0.0), xpr(0.0), xpp(0.0) {};
    EnergyDens(RDFType rtype, RDF& rr, RDF& rp, RDF& pr, RDF& pp) {
      xrr=rr.pint(rtype);
      xrp=rp.pint(rtype);
      xpr=pr.pint(rtype); 
      xpp=pp.pint(rtype);
    }
    void calc(double rdens, double pdens, Parameter& p, double& murna, double& mupos);
    void calc(double rdens, double pdens, Parameter& p, double vol, double& murna, double& mupos);
    void calc(double rdens, double pdens, double erdens, double epdens, Parameter& p, double& murna, double& mupos);
};

double EnergyDens::svol(double vratio, double temp) {
  if (vratio<=0.0000000001) {
    return 0.0;
  } else {
    return -temp*kb*log(vratio);
  }
}

void EnergyDens::calc(double rdens, double pdens, Parameter& p, double vol, double& murna, double& mupos) {
  double hrna=(xrr*rdens+xrp*pdens)/2.0;
  double hpos=(xpp*pdens+xpr*rdens)/2.0;
  double srna=svol(vol/p.sysvol,p.temp);
  double spos=svol(vol/p.sysvol,p.temp);
  murna=hrna+srna;
  mupos=hpos+spos;
} 

void EnergyDens::calc(double rdens, double pdens, Parameter& p, double& murna, double& mupos) {
  double hrna=(xrr*rdens+xrp*pdens)/2.0;
  double hpos=(xpp*pdens+xpr*rdens)/2.0;
  double srna=svol(p.rnadens()/rdens,p.temp);
  double spos=svol(p.posdens()/pdens,p.temp);
  murna=hrna+srna;
  mupos=hpos+spos;
} 

void EnergyDens::calc(double rdens, double pdens, double erdens, double epdens, Parameter& p, double& murna, double& mupos) {
  double hrna=(xrr*rdens+xrp*pdens)/2.0;
  double hpos=(xpp*pdens+xpr*rdens)/2.0;
  double srna=svol(p.rnadens()/erdens,p.temp);
  double spos=svol(p.posdens()/epdens,p.temp);
  murna=hrna+srna;
  mupos=hpos+spos;
} 

double gmix(double n, double x, double t) {
  return n*kb*x*log(x)*t;
}

double radius(double vol) {
  if (vol<=0.0) 
     return 0.0;
  else 
     return pow(3.0*vol/4.0/pi,1.0/3.0);
}

double concentration(double density) {
  return density/navo/1.0E-24*1000.0; // result in mM
}

double numberdensity(double conc) {
  return conc/1000.0*1.0E-24*navo;
}

void advance(double &n, double scale=0.02) {
  double ln=log(n)/l10;
  ln+=scale;
  n=pow(10,ln);
}

void mucalc(double nr, double np, double clustvol, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster,
            double &rld, double &rhd, double &pld, double &phd, double &murna, double &mupos, double &dr, double &dp, double &g) {

  double effcvol=clustvol-(par.nrna-nr)*par.volrna-(par.npos-np)*par.volpos;

  rld=nr/(par.sysvol-clustvol);
  rhd=((double)par.nrna-nr)/clustvol;
  double rhde=((double)par.nrna-nr)/effcvol;

  pld=np/(par.sysvol-clustvol);
  phd=((double)par.npos-np)/clustvol;
  double phde=((double)par.npos-np)/effcvol;

  double murnalow=0.0;
  double murnahi=0.0;
  double muposlow=0.0;
  double muposhi=0.0;

  edilute.calc(rld,pld,par,murnalow,muposlow);
  ecluster.calc(rhd,phd,rhde,phde,par,murnahi,muposhi);

  dr=fabs(murnahi-murnalow);
  dp=fabs(muposhi-muposlow);

  murna=(murnalow+murnahi)/2.0;
  mupos=(muposlow+muposhi)/2.0;

  g=murna*par.nrna+mupos*par.npos;
}


double npfunc(double nr, double np, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double xdrr=edilute.xrr;
  double xdrp=edilute.xrp;
  double xcrr=ecluster.xrr;
  double xcrp=ecluster.xrp;

  double v=par.sysvol;
  double vl=v-vc;
  double nrc=par.nrna-nr;
  double npc=par.npos-np;

  double effcvol=vc-nrc*par.volrna-npc*par.volpos;
  
  double val=xdrr*nr/vl-xcrr*nrc/vc-xcrp*par.npos/vc;
  val+=np*(xdrp/vl+xcrp/vc);
  val-=par.kt2*log(vl*nrc/nr/effcvol);

  return val;
}

double nrfunc(double nr, double np, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double xdpr=edilute.xpr;
  double xdpp=edilute.xpp;
  double xcpr=ecluster.xpr;
  double xcpp=ecluster.xpp;

  double v=par.sysvol;
  double vl=v-vc;
  double nrc=par.nrna-nr;
  double npc=par.npos-np;

  double effcvol=vc-nrc*par.volrna-npc*par.volpos;

  double val=xdpp*np/vl-xcpp*npc/vc-xcpr*par.nrna/vc;
  val+=nr*(xdpr/vl+xcpr/vc);
  val-=par.kt2*log(vl*npc/np/effcvol);

  return val;
}

double npderiv(double nr, double np, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double xdrp=edilute.xrp;
  double xcrp=ecluster.xrp;

  double v=par.sysvol;

  double effcvol=vc-(par.nrna-nr)*par.volrna-(par.npos-np)*par.volpos;

  double val=xdrp/(v-vc)+xcrp/vc+par.kt2*par.volpos/effcvol; 

  return val;
}

double nrderiv(double nr, double np, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double xdpr=edilute.xpr;
  double xcpr=ecluster.xpr;

  double v=par.sysvol;

  double effcvol=vc-(par.nrna-nr)*par.volrna-(par.npos-np)*par.volpos;

  double val=xdpr/(v-vc)+xcpr/vc+par.kt2*par.volrna/effcvol; 

  return val;
}

double npderiv0(double nr, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double xdrp=edilute.xrp;
  double xcrp=ecluster.xrp;

  double v=par.sysvol;

  double val=-par.kt2/(xdrp/(v-vc)+xcrp/vc)-vc/par.volpos+(par.nrna-nr)*par.volrna/par.volpos+par.npos;
  return val;
} 

double nrderiv0(double np, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double xdpr=edilute.xpr;
  double xcpr=ecluster.xpr;

  double v=par.sysvol;

  double val=-par.kt2/(xdpr/(v-vc)+xcpr/vc)-vc/par.volrna+(par.npos-np)*par.volpos/par.volrna+par.nrna;
  return val;
} 

double nrdnpzero(double nr, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double xdrr=edilute.xrr;
  double xcrr=ecluster.xrr;

  double v=par.sysvol;

  double effcvol=vc-(par.nrna-nr)*par.volrna-(par.npos)*par.volpos;
  double val=xdrr/(v-vc)+xcrr/vc+par.kt2*(par.nrna/((par.nrna-nr)*nr)+par.volrna/effcvol);
  return val;
} 

double npdnrzero(double np, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double xdpp=edilute.xpp;
  double xcpp=ecluster.xpp;

  double v=par.sysvol;

  double effcvol=vc-(par.nrna)*par.volrna-(par.npos-np)*par.volpos;
  double val=xdpp/(v-vc)+xcpp/vc+par.kt2*(par.npos/((par.npos-np)*np)+par.volpos/effcvol);
  return val;
}
 
double newraphnp(double inp, double min, double max, double nr, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double np=inp;
  double v=npfunc(nr,np,vc,par,edilute,ecluster);
  int nrun=0;
  do {
    double dv=npderiv(nr,np,vc,par,edilute,ecluster);
    np=(-v+np*dv)/dv;
    if (np<min) np=min;
    if (np>max) np=max;
    v=npfunc(nr,np,vc,par,edilute,ecluster);
  } while (fabs(v)>0.000000001 && ++nrun<100);

  if (fabs(v)>0.00001) np=-1.0;

  return np;
}

double newraphnr(double inr, double min, double max, double np, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double nr=inr;
  double v=nrfunc(nr,np,vc,par,edilute,ecluster);
  int nrun=0;
  do {
    double dv=nrderiv(nr,np,vc,par,edilute,ecluster);
    nr=(-v+nr*dv)/dv;
    if (nr<min) nr=min;
    if (nr>max) nr=max;
    v=nrfunc(nr,np,vc,par,edilute,ecluster);
  } while (fabs(v)>0.000000001 && ++nrun<100);

  if (fabs(v)>0.00001) nr=-1.0;

  return nr;
}

double newraphnrnpzero(double inr, double min, double max, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double nr=inr;
  double v=npfunc(nr,0.0,vc,par,edilute,ecluster);
  int nrun=0;
  do {
    double dv=nrdnpzero(nr,vc,par,edilute,ecluster);
    nr=(-v+nr*dv)/dv;
    if (nr<min) nr=min;
    if (nr>max) nr=max;
    v=npfunc(nr,0.0,vc,par,edilute,ecluster);
  } while (fabs(v)>0.000000001 && ++nrun<100);

  if (fabs(v)>0.00001) nr=-1.0;

  return nr;
}

double newraphnpnrzero(double inp, double min, double max, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double np=inp;
  double v=nrfunc(0.0,np,vc,par,edilute,ecluster);
  int nrun=0;
  do {
    double dv=npdnrzero(np,vc,par,edilute,ecluster);
    np=(-v+np*dv)/dv;
    if (np<min) np=min;
    if (np>max) np=max;
    v=nrfunc(0.0,np,vc,par,edilute,ecluster);
  } while (fabs(v)>0.000000001 && ++nrun<100);
  
  if (fabs(v)>0.00001) np=-1.0;

  return np;
}

double solvenp(double nr, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double minnp=(vc-0.20*vc-(par.nrna-nr)*par.volrna-par.npos*par.volpos)/(-par.volpos);

  if (minnp>par.npos) { return -1; }
  if (minnp<0) { minnp=0.0; }

  double minmax=npderiv0(nr,vc,par,edilute,ecluster);

  double np=-1.0;
  if (minmax<minnp || minmax>par.npos) {
    double v0=npfunc(nr,minnp,vc,par,edilute,ecluster);
    double vmax=npfunc(nr,par.npos,vc,par,edilute,ecluster);
    if ((v0>0 && vmax<0) || (v0<0 && vmax>0)) {          // straight downhill or uphill
       np=newraphnp(minnp,minnp,par.npos,nr,vc,par,edilute,ecluster);
    }
  } else if (minmax>minnp && minmax<par.npos) {          // maximum/minimum in target range
    double v0=npfunc(nr,minnp,vc,par,edilute,ecluster);
    double vminmax=npfunc(nr,minmax,vc,par,edilute,ecluster);
    double vmax=npfunc(nr,par.npos,vc,par,edilute,ecluster);

    double np1=-1.0;
    double np2=-1.0;
    if ((v0<0 && vminmax>=0) || (v0>0 && vminmax<=0)) {
      np1=newraphnp(minnp,minnp,par.npos,nr,vc,par,edilute,ecluster);
    } 
    if ((vminmax>0 && vmax<0) || (vminmax<0 && vmax>0)) {
      np2=newraphnp(par.npos,minnp,par.npos,nr,vc,par,edilute,ecluster);
    }
    if (np1>0 && np1<par.npos) {
      np=np1;
      if (np2>0 && np2<par.npos) {
        if (verbose) fprintf(stderr,"two solutions: %lf %lf\n",np1,np2);
      }
    } else if (np2>0 && np2<par.npos) {
      np=np2;
    } 
  }

  return np;
}

double solvenr(double np, double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double minnr=(vc-0.20*vc-(par.npos-np)*par.volpos-par.nrna*par.volrna)/(-par.volrna);

  if (minnr>par.nrna) { return -1; }
  if (minnr<0) { minnr=0.0; }

  double minmax=nrderiv0(np,vc,par,edilute,ecluster);

  double nr=-1.0;
  if (minmax<minnr || minmax>par.nrna) {
    double v0=nrfunc(minnr,np,vc,par,edilute,ecluster);
    double vmax=nrfunc(par.nrna,np,vc,par,edilute,ecluster);
    if ((v0>0 && vmax<0) || (v0<0 && vmax>0)) {       // straight downhill or uphill
       nr=newraphnr(minnr,minnr,par.nrna,np,vc,par,edilute,ecluster);
    }
  } else if (minmax>minnr && minmax<par.nrna) {       // maximum/minimum in target range
    double v0=nrfunc(minnr,np,vc,par,edilute,ecluster);
    double vminmax=nrfunc(minmax,np,vc,par,edilute,ecluster);
    double vmax=nrfunc(par.nrna,np,vc,par,edilute,ecluster);

    double nr1=-1.0;
    double nr2=-1.0;
    if ((v0<0 && vminmax>=0) || (v0>0 && vminmax<=0)) {
      nr1=newraphnr(minnr,minnr,par.nrna,np,vc,par,edilute,ecluster);
    } 
    if ((vminmax>0 && vmax<0) || (vminmax<0 && vmax>0)) {
      nr2=newraphnr(par.nrna,minnr,par.nrna,np,vc,par,edilute,ecluster);
    }
    if (nr1>0 && nr1<par.nrna) {
      nr=nr1;
      if (nr2>0 && nr2<par.nrna) {
        fprintf(stderr,"two solutions: %lf %lf\n",nr1,nr2);
      }
    } else if (nr2>0 && nr2<par.nrna) {
      nr=nr2;
    } 
  }

  return nr;
}

double solvenrnpzero(double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double minnr=(vc-0.20*vc-(par.npos)*par.volpos-par.nrna*par.volrna)/(-par.volrna);

  if (minnr>par.nrna) { return -1; }
  if (minnr<0.001) { minnr=0.001; }

  double nr=-1.0;
  double v0=npfunc(minnr,0.0,vc,par,edilute,ecluster);
  double vmax=npfunc(par.nrna*0.99,0.0,vc,par,edilute,ecluster);

  if ((v0>0 && vmax<0) || (v0<0 && vmax>0)) {       // straight downhill or uphill
     nr=newraphnrnpzero(minnr,minnr,par.nrna,vc,par,edilute,ecluster);
  }
  return nr;
}

double solvenpnrzero(double vc, Parameter &par, EnergyDens &edilute, EnergyDens &ecluster) {
  double minnp=(vc-0.20*vc-(par.nrna)*par.volrna-par.npos*par.volpos)/(-par.volpos);

  if (minnp>par.npos) { return -1; }
  if (minnp<0.000001) { minnp=0.000001; }

  double np=-1.0;
  double v0=nrfunc(0.0,minnp,vc,par,edilute,ecluster);
  double vmax=nrfunc(0.0,par.npos*0.999,vc,par,edilute,ecluster);

  if ((v0>0 && vmax<0) || (v0<0 && vmax>0)) {       // straight downhill or uphill
     np=newraphnpnrzero(minnp,minnp,par.npos,vc,par,edilute,ecluster);
  }
  return np;
}

int main(int argc, char **argv) {
  Parameter par;

  double rdfcut=20.0;    // g(r) beyond this cutoff is set to 1.0
  double maxrad=100.0;   // maximum radius for integration 
  double thresh=0.0005;  // threshold to determine two energies are the same

  char tag[512];
  strcpy(tag,"");

  char rnarnardfname[512];
  char rnaposrdfname[512];
  char posrnardfname[512];
  char posposrdfname[512];
  
  strcpy(rnarnardfname,"trna.trna"); 
  strcpy(rnaposrdfname,"trna.pos2"); 
  strcpy(posrnardfname,"pos2.trna"); 
  strcpy(posposrdfname,"pos2.pos2"); 

  for (int i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-verbose")) {
      verbose=1;
    } else if (!strcmp(argv[i],"-highres")) {
      highres=1;
    } else if (!strcmp(argv[i],"-rnarna")) {
      strcpy(rnarnardfname,argv[++i]);
    } else if (!strcmp(argv[i],"-rnapos")) {
      strcpy(rnaposrdfname,argv[++i]);
    } else if (!strcmp(argv[i],"-posrna")) {
      strcpy(posrnardfname,argv[++i]);
    } else if (!strcmp(argv[i],"-pospos")) {
      strcpy(posposrdfname,argv[++i]);
    } else if (!strcmp(argv[i],"-tag")) {
      strcpy(tag,argv[++i]);
    } else if (!strcmp(argv[i],"-rdfcut")) {
      rdfcut=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-maxrad")) {
      maxrad=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-thresh")) {
      thresh=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-temp")) {
      par.settemp(atof(argv[++i]));
    } else if (!strcmp(argv[i],"-q")) {
      Field f(argv[++i],(char *)":");
      par.qrna=atof(f[0]);
      par.qpos=atof(f[1]);
    } else if (!strcmp(argv[i],"-qrna")) {
      par.qrna=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-qpos")) {
      par.qpos=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-r")) {
      Field f(argv[++i],(char *)":");
      par.setrrna(atof(f[0]));
      par.setrpos(atof(f[1]));
    } else if (!strcmp(argv[i],"-rrna")) {
      par.setrrna(atof(argv[++i]));
    } else if (!strcmp(argv[i],"-rpos")) {
      par.setrpos(atof(argv[++i]));
    } else if (!strcmp(argv[i],"-kappa")) {
      par.kappa=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-vfac")) {
      par.vfac=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-nrscale")) {
      par.nrscale=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-a0")) {
      par.a0=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-epsilon")) {
      par.epsilon=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-len")) {
      par.setsyslen(atof(argv[++i]));
    } else if (!strcmp(argv[i],"-n")) {
      Field f(argv[++i],(char *)":");
      par.nrna=atof(f[0]);
      par.npos=atof(f[1]);
    } else if (!strcmp(argv[i],"-nrna")) {
      par.nrna=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-npos")) {
      par.npos=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-c")) {
      Field f(argv[++i],(char *)":");
      double crna=atof(f[0]);
      double cpos=atof(f[1]);
      par.nrna=numberdensity(crna)*par.sysvol;
      par.npos=numberdensity(cpos)*par.sysvol;
    } else if (!strcmp(argv[i],"-crna")) {
      double crna=atof(argv[++i]);
      par.nrna=numberdensity(crna)*par.sysvol;
    } else if (!strcmp(argv[i],"-cpos")) {
      double cpos=atof(argv[++i]);
      par.npos=numberdensity(cpos)*par.sysvol;
    }   
  }

  RDF rnarna(rnarnardfname,par.sigrna(),par.dsigrna(),par.arna()*par.arna(),par,rdfcut,maxrad);
  RDF rnapos(rnaposrdfname,(par.sigrna()+par.sigpos())/2.0,(par.dsigrna()+par.dsigpos())/2.0,
                           par.arna()*par.apos(),par,rdfcut,maxrad);
  RDF posrna(posrnardfname,(par.sigpos()+par.sigrna())/2.0,(par.dsigpos()+par.dsigrna())/2.0,
                           par.apos()*par.arna(),par,rdfcut,maxrad);
  RDF pospos(posposrdfname,par.sigpos(),par.dsigpos(),par.apos()*par.apos(),par,rdfcut,maxrad);

  EnergyDens edisperse(DISPERSE,rnarna,rnapos,posrna,pospos);
  EnergyDens edilute(DILUTE,rnarna,rnapos,posrna,pospos);
  EnergyDens ecluster(CLUSTER,rnarna,rnapos,posrna,pospos);

  double murnadisperse,muposdisperse;
  edisperse.calc(par.rnadens(),par.posdens(),par,murnadisperse,muposdisperse);
  double gdisperse=murnadisperse*par.nrna+muposdisperse*par.npos;
 
  double gmixall=gmix(par.nrna+par.npos,par.nrna/(par.nrna+par.npos),par.temp)+
                 gmix(par.nrna+par.npos,par.npos/(par.nrna+par.npos),par.temp); 

  gdisperse+=gmixall;

  double ntot=(double)par.nrna+(double)par.npos;

  double bestnr=-1;
  double bestnp=-1;
  double bestcvol=-1;
  double ming=9999999;

  double minvolclusterall=-1;
  double gminclusterall=999999;

  double minvolclusternpzero=-1;
  double gminclusternpzero=999999;
  double minnrnpzero=-1;

  double minvolclusternrzero=-1;
  double gminclusternrzero=999999;
  double minnpnrzero=-1;

  double clustvol=0.3*par.sysvol;
  do {
    double murnaclusterall,muposclusterall;
    double effvol=clustvol-par.nrna*par.volrna-par.npos*par.volpos;
    ecluster.calc(par.nrna/clustvol,par.npos/clustvol,par.nrna/effvol,par.npos/effvol,par,murnaclusterall,muposclusterall);
    double gclusterall=murnaclusterall*par.nrna+muposclusterall*par.npos;
    gclusterall+=gmixall;
    if (gclusterall<gminclusterall && (par.volrna/(clustvol/par.nrna)+par.volpos/(clustvol/par.npos))<0.5) {
       gminclusterall=gclusterall;
       minvolclusterall=clustvol;
    }

    double nrnpzero=solvenrnpzero(clustvol,par,edilute,ecluster);
    if (nrnpzero>0 && nrnpzero<par.nrna) {
      effvol=clustvol-(par.nrna-nrnpzero)*par.volrna-par.npos*par.volpos;
      double murnaclusternpzero,muposclusternpzero;
      ecluster.calc((par.nrna-nrnpzero)/clustvol,par.npos/clustvol,(par.nrna-nrnpzero)/effvol,par.npos/effvol,par,murnaclusternpzero,muposclusternpzero);
      double murnadilutenpzero,muposdilutenpzero;
      edilute.calc(nrnpzero/(par.sysvol-clustvol),0.0,par,murnadilutenpzero,muposdilutenpzero);
      double gclusternpzero=murnaclusternpzero*par.nrna+muposclusternpzero*par.npos;

      double gmixt=gmix((par.nrna-nrnpzero)+par.npos,par.npos/((par.nrna-nrnpzero)+par.npos),par.temp)+
                   gmix((par.nrna-nrnpzero)+par.npos,(par.nrna-nrnpzero)/((par.nrna-nrnpzero)+par.npos),par.temp); 

      gclusternpzero+=gmixt;

      if (gclusternpzero<gminclusternpzero && (par.volrna/(clustvol/(par.nrna-nrnpzero))+par.volpos/(clustvol/par.npos))<0.5) {
         gminclusternpzero=gclusternpzero;
         minvolclusternpzero=clustvol;
         minnrnpzero=nrnpzero;
      }
    }
     
    double npnrzero=solvenpnrzero(clustvol,par,edilute,ecluster);
    if (npnrzero>0 && npnrzero<par.npos) {
      effvol=clustvol-(par.nrna)*par.volrna-(par.npos-npnrzero)*par.volpos;
      double murnaclusternrzero,muposclusternrzero;
      ecluster.calc(par.nrna/clustvol,(par.npos-npnrzero)/clustvol,par.nrna/effvol,(par.npos-npnrzero)/effvol,par,murnaclusternrzero,muposclusternrzero);
      double murnadilutenrzero,muposdilutenrzero;
      edilute.calc(0.0,npnrzero/(par.sysvol-clustvol),par,murnadilutenrzero,muposdilutenrzero);
      double gclusternrzero=murnaclusternrzero*par.nrna+muposclusternrzero*par.npos;

      double gmixt=gmix(par.nrna+(par.npos-npnrzero),par.nrna/(par.nrna+(par.npos-npnrzero)),par.temp)+
                   gmix(par.nrna+(par.npos-npnrzero),(par.npos-npnrzero)/(par.nrna+(par.npos-npnrzero)),par.temp); 

      gclusternrzero+=gmixt;

      if (gclusternrzero<gminclusternrzero && (par.volrna/(clustvol/par.nrna)+par.volpos/(clustvol/(par.npos-npnrzero)))<0.5) {
         gminclusternrzero=gclusternrzero;
         minvolclusternrzero=clustvol;
         minnpnrzero=npnrzero;
      }
    }

    double initialnr=(double)par.nrna/100000.0;
    double nr=initialnr;
    int found=0;
    do {
        double np=solvenp(nr,clustvol,par,edilute,ecluster);

        if (np>0 && (np<par.npos/2 || nr<par.nrna/2)) {
          double rld,rhd,pld,phd,dr,dp,g,murna,mupos;
          mucalc(nr,np,clustvol,par,edilute,ecluster,rld,rhd,pld,phd,murna,mupos,dr,dp,g);  

          if (dr<thresh && dp<thresh) {
            found=1;

            double gnr=solvenr(np,clustvol,par,edilute,ecluster);

            double grld,grhd,gpld,gphd,gdr,gdp,gg,gmurna,gmupos;
            mucalc(gnr,np,clustvol,par,edilute,ecluster,grld,grhd,gpld,gphd,gmurna,gmupos,gdr,gdp,gg);  
            
            if (gdr*gdr+gdp*gdp<dr*dr+dp*dp) {
              nr=gnr;
              rld=grld;
              rhd=grhd;
              pld=gpld;
              phd=gphd;
              murna=gmurna;
              mupos=gmupos;
              dr=gdr;
              dp=gdp;
              g=gg;
            }
 
            double gmixt=gmix(nr+np,nr/(nr+np),par.temp)+gmix(nr+np,np/(nr+np),par.temp)+
                         gmix((par.nrna-nr)+(par.npos-np),(par.nrna-nr)/((par.nrna-nr)+(par.npos-np)),par.temp)+
                         gmix((par.nrna-nr)+(par.npos-np),(par.npos-np)/((par.nrna-nr)+(par.npos-np)),par.temp);

            g+=gmixt;

            if (verbose) {
                fprintf(stderr,"%1.10lf (%1.2lf)  %1.10lf (%1.2lf) | %1.10lf (%1.2lf) %1.10lf (%1.2lf) :: %lf %lf :: %lf %lf :: %lf %lf :: %lf %lf ::: %lf\n",
                              rld,nr,rhd,(par.nrna-nr),pld,np,phd,par.npos-np,
                              murna,mupos,dr,dp,gdisperse,g,par.volrna/(1.0/rhd),par.volpos/(1.0/phd),clustvol);
            }

            if (g<ming && (par.volrna/(1.0/rhd)+par.volpos/(1.0/phd))<0.5) {
              ming=g;
              bestnr=nr;
              bestnp=np; 
              bestcvol=clustvol;
            }
          }
        }
        advance(nr,par.nrscale);
     } while (nr<par.nrna && !found);   
     clustvol*=par.vfac;
  } while (clustvol>=0.001*par.sysvol); 

  char solution[100]; 
  double rld,rhd,pld,phd,cvol;
  if (gdisperse<=ming && gdisperse<=gminclusternrzero && gdisperse<=gminclusternpzero && gdisperse<=gminclusterall) {
    rld=rhd=par.rnadens();
    pld=phd=par.posdens();
    cvol=par.sysvol;
    strcpy(solution,"disperse");
  } else if (ming<=gminclusternrzero && ming<=gminclusternpzero && ming<=gminclusterall && ming<=gdisperse) {
    cvol=bestcvol;
    rld=bestnr/(par.sysvol-cvol);
    rhd=((double)par.nrna-bestnr)/cvol;
    pld=bestnp/(par.sysvol-cvol);
    phd=((double)par.npos-bestnp)/cvol;
    strcpy(solution,"phasesep");
  } else if (gminclusternrzero<=ming && gminclusternrzero<=gminclusternpzero && gminclusternrzero<=gminclusterall && gminclusternrzero<=gdisperse) {
    cvol=minvolclusternrzero;
    rld=0.0;
    rhd=((double)par.nrna)/cvol;
    pld=minnpnrzero/(par.sysvol-cvol);
    phd=((double)par.npos-minnpnrzero)/cvol;
    strcpy(solution,"phasesep[nr0]");
  } else if (gminclusternpzero<=ming && gminclusternpzero<=gminclusternrzero && gminclusternpzero<=gminclusterall && gminclusternpzero<=gdisperse) {
    cvol=minvolclusternpzero;
    rld=minnrnpzero/(par.sysvol-cvol);
    rhd=((double)par.nrna-minnrnpzero)/cvol;
    pld=0.0;
    phd=((double)par.npos)/cvol;
    strcpy(solution,"phasesep[np0]");
  } else if (gminclusterall<=gminclusternrzero && gminclusterall<=gminclusternpzero && gminclusterall<=ming && gminclusterall<=gdisperse) {
    cvol=minvolclusterall;
    rld=0.0;
    rhd=((double)par.nrna)/cvol;
    pld=0.0;
    phd=((double)par.npos)/cvol;
    strcpy(solution,"condensed");
  } else {
    cvol=100;
    rld=0.0;
    rhd=0.0;
    pld=0.0;
    phd=0.0;
    strcpy(solution,"none");
  }

  if (highres) {
     fprintf(stdout,"%s %s :[rna_mM]: %1.12lf %1.12lf :[pos_mM]: %1.12lf %1.12lf :[clusterrad_nm]: %1.12lf :[volfrac]: %1.12lf %1.12lf\n",tag,solution,
             concentration(rld),concentration(rhd),concentration(pld),concentration(phd),radius(cvol),par.volrna/(1.0/rhd),par.volpos/(1.0/phd));
  } else {
     fprintf(stdout,"%-5s %-14s :[rna_mM]: %7.6lf %5.2lf :[pos_mM]: %7.6lf %5.2lf :[clusterrad_nm]: %5.2lf :[volfrac]: %1.2lf %1.2lf\n",tag,solution,
             concentration(rld),concentration(rhd),concentration(pld),concentration(phd),radius(cvol),par.volrna/(1.0/rhd),par.volpos/(1.0/phd));
  }
}     


/*
 * var.cpp
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */
#include "var.h"


namespace variable {
  
  //******************************************************************************** Operators
  
  var operator*(const double& scalar,const var& var1){ // Pre-multiplication with scalar
    TRACE(0,"operator*(scalar,var)");
    assert(var1.gc!=NULL);
    vd newtdata=scalar*var1.tdata();
    return var(*(var1.gc),newtdata);
  }

  //***************************************** The var class
  var::var(const Globalconf& gc): var(gc,0.0) {  }
  var::var(const Globalconf& gc,double initval) :gc(&gc),Nf(gc.Nf),Ns(gc.Ns) {
    TRACE(0,"var::var(const Globalconf& gc, double initval)");
    timedata=vd(Ns);
    amplitudedata=vd(Ns);
    settdata(initval);
  }
  var::var(const Globalconf& gc,const vd& timedata):var(gc){ // Create a variable and fill it with time data.
    TRACE(0,"var::var(gc,timedata)");
    this->timedata=timedata;
    dft();
  }
  void var::updateNf(){
    TRACE(0,"var::updateNf()");
    if(this->Ns!=gc->Ns){
      TRACE(5,"UPDATENF untested code");
      assert((Ns%2)==1);	// Check if number of samples is not even
      if(this->Ns>gc->Ns){
	amplitudedata=amplitudedata.subvec(0,gc->Ns-1);
      }
      if(this->Ns<gc->Ns){
	vd oldadata=amplitudedata;
	amplitudedata=vd(gc->Ns,fillwith::zeros);
	amplitudedata.subvec(0,this->Ns-1)=oldadata;
      }
      this->Ns=gc->Ns;		       // Update this number of samples
      timedata=vd(Ns,fillwith::zeros); // Reinitialize timedata
      idft();
    }
  }
  // var var::operator()(const var& v) { //Copy constructor
  //   TRACE(0,"var copy constructor");
  //   assert(v.gc!=NULL);
  //   vd tdata=v.tdata();
  //   return var(*(this->gc),tdata);
  // }
  // var& var::operator=(const var& v){
  //   TRACE(0,"var::operator=(const var& v)");
  //   this->timedata=v.timedata;
  //   //	TRACE(0,"this->timedata:" <<this->timedata);
  //   //	TRACE(0,"v.timedata:" <<v.timedata);
  //   this->amplitudedata=v.amplitudedata;
  //   return *this;
  // }
  var var::operator*(const var& var2) const { // Multiply two variables in timed
    TRACE(0,"var::operator*(const var& var2) const");
    assert(this->Ns==var2.Ns);
    vd tdata=this->tdata()%var2.tdata();
    return var(*(this->gc),tdata);
  }
  var var::operator*(const d& scalar) const {	// Post-multiplication with scalar
    assert(this->gc!=NULL);
    vd thisadata=this->tdata();
    return var(*(this->gc),scalar*thisadata);
  }
  // Get methods (which require implementation)
  const d& var::operator()(us i) const {//Extract result at specific frequency
      TRACE(0,"var::operator(us i)");
      assert(i<Ns);
      TRACE(-1,"amplitudedata: "<<amplitudedata);
      return amplitudedata(i);
  }
  d& var::operator()(us i) {//Extract result at specific frequency
      TRACE(0,"var::operator(us i)");
      assert(i<Ns);
      TRACE(-1,"amplitudedata: "<<amplitudedata);
      return amplitudedata(i);
    }  
  vc var::getcRes() const
  {
    TRACE(0,"var::getcRes()");
    //	TRACE(0,"amplitudedata:" << amplitudedata);
    vc cadata(Nf+1);
    cadata(0)=amplitudedata(0);
    for(us i=1;i<Nf+1;i++)
      cadata(i)=amplitudedata(2*i-1)-I*amplitudedata(2*i); //The minus is very important
    //	TRACE(0,"Resulting cadata:" << cadata);
    return cadata;
  }
  d var::tdata(d t) const //Extract the value for a given time
  {
    d result=0;
    vc cres=getcRes();
    for(us n=0;n<Nf+1;n++)
      {
	result+=real(cres(n)*exp(I*double(n)*gc->omg*t));
      }
    return result;
  }

  // Set methods
  void var::set(double val,us freqnr) { //Set result for specific frequency zero,real one, -imag one, etc
    updateNf();
    amplitudedata[freqnr]=val;
    idft();
    TRACE(-3,"var::set(d val,us freqnr) adata:"<<amplitudedata);
  }
  void var::set(const vc res)
  {
    TRACE(0,"var::set(const vc& res)");
    updateNf();
    amplitudedata(0)=res(0).real();
    for(us i=1;i<Nf+1;i++){
      amplitudedata(2*i-1)=res(i).real();
      amplitudedata(2*i)=-1.0*res(i).imag();
    }
    idft();
  }
  void var::set(const vd val) {
    TRACE(0,"var::set(const vd& val)");
    updateNf();
    amplitudedata=val;
    idft();
  }
  void var::setResfluc(vd& val) {
    amplitudedata.subvec(1,Ns-1)=val;
  }
  void var::settdata(double val) {
    TRACE(0,"var::settdata(double val)");
    updateNf();
    timedata.fill(val);
    dft();
  }
  void var::settdata(vd& val) {
    TRACE(0,"var::settdata(vd& val)");
    updateNf();
    assert(val.size()==timedata.size());
    timedata=val;
    dft();
  }
  //Show methods
  void var::showtdata() const {
    unsigned i;
    cout << "[" ;
    for(i=0; i<Ns-1; i++) {
      cout << timedata[i] << " ";
    }
    cout << timedata[Ns-1] << "]\n";

  }
  void var::showRes() const {
    unsigned i;
    cout << "[" ;
    for(i=0; i<Ns-1; i++) {
      cout << amplitudedata[i] << " ";
    }
    cout << amplitudedata[Ns-1] << "]\n";

  }

  // Internal methods for syncing time and amplitude data
  void var::dft() {
    TRACE(0,"var::dft()");
    amplitudedata=gc->fDFT*timedata;
  }
  void var::idft() { //Internal idft
    timedata=gc->iDFT*amplitudedata;
  }

  //Get a variable which is the time derivative of the current one
  var var::ddt() const {
    var result(*(this->gc));
    vd newadata=gc->DDTfd*amplitudedata;
    result.set(newadata);
    return result;
  }
  // The product

  var var::operator/(const var& var2) const
  {
    vd tdata=this->tdata()/var2.tdata();
    var newvar(*(this->gc));
    newvar.settdata(tdata);
    return newvar;
  }
  var::~var() {// The destructor
    TRACE(-5,"var destructor called");
  }
  //***************************************** End of the var class


  ostream& operator<< (ostream& out,var& v){
    vd res=v();
    return out << res;
  }


  //The vvar vlass
  // vvar::vvar() {	}
  // vvar::vvar(string name,us gp,us Nf): vvar(name,gp,Nf,0.0)  {}
  // vvar::vvar(string name,us gp,us Nf,double initv) :
  //   Name(name),gp(gp),Nf(Nf) {
  //   Ns=2*Nf+1; //Number of time samples
  //   Dofs=Ns*gp;
  //   TRACE(0,"Warning, vvar not changed for Globalconf!!");
  //   //for (us i=0;i<gp;i++) data.push_back(var(Nf,initv));
  //   Error=vd(Dofs);
  //   Result=vd(Dofs);
  // }
  // string vvar::name() {return Name;}
  // vd& vvar::getRes() {
  //   //Create a column with all data ordered
  //   for(us i=0;i<gp;i++) {
  //     vd res=data[i]();
  //     //cout << t.size() << endl;
  //     //cout << res.size() << endl;
  //     Result.subvec(i*Ns,(i+1)*Ns-1)=res;
  //   }
  //   return Result;
  // }
  // vd vvar::getRes(us freq) {
  //   vd Res=vd(gp);
  //   for(us i=0;i<gp;i++) {
  //     Res[i]=data[i](freq);
  //   }
  //   return Res;
  // }
  // void vvar::setRes(vd& Res) {
  //   for(us i=0;i<gp;i++){
  //     vd resi=Res.subvec(i*Ns,(i+1)*Ns-1);
  //     data[i].set(resi);
  //   }
  // }
  // void vvar::setRes(vd& Res,us freqnr) {
  //   for(us i=0;i<gp;i++){
  //     data[i].set(Res[i],freqnr);
  //   }
  // }
  // void vvar::showResult() {
  //   unsigned i;
  //   vd& result=getRes();
  //   cout << "Total data for " << Name << endl;
  //   cout << "[ " ;
  //   for (i=0;i<Dofs-1;i++){
  //     cout << result[i] << " ";
  //   }
  //   cout << result[Dofs-1] << " ]\n";
  // }


  // //	vd vvar::dotn(us i,vd x) {
  // //		// dot n. For the boundaries, extrapolation is used
  // //		vd vi = data[i].getAdata();
  // //		vd vi_plus_half;
  // //		vd vi_min_half;
  // //		vd vip1;
  // //		vd vim1;
  // //		if(i>0 &&i<gp-1) {
  // //			vip1=data[i+1].getAdata();
  // //			vim1=data[i-1].getAdata();
  // //			vi_plus_half=0.5*(vip1+vi);
  // //			vi_min_half=0.5*(vi+vim1);
  // //		}
  // //		else if(i==0){
  // //			vip1=data[i+1].getAdata();
  // //			vi_plus_half=0.5*(vip1+vi);
  // //			vd dvdx_plus_half=(vip1-vi)/(x[i+1]-x[i]);
  // //			double h=(x[i+1]-x[i]);
  // //			vi_min_half=vi_plus_half-dvdx_plus_half*h;
  // //		}
  // //		else {
  // //			vim1=data[i-1].getAdata();
  // //			vi_min_half=0.5*(vim1+vi);
  // //			double h=(x[i]-x[i-1]);
  // //			vd dvdx_i_min_half=(vi-vim1)/h;
  // //			vi_plus_half=vi_min_half+dvdx_i_min_half*h;
  // //		}
  // //	return vi_plus_half-vi_min_half;
  // //	}
  // vd vvar::ddx_forward(us i,const vd& x){
  //   vd result(Ns);
  //   assert((i>=0) && (i<=gp));
  //   if (i<gp-1){
  //     result=(data[i+1]()-data[i]())/(x[i+1]-x[i]);
  //     return result;
  //   }
  //   else{
  //     result=(data[i]()-data[i-1]())/(x[i]-x[i-1]);
  //     return result;
  //   }

  // }
  // vd vvar::ddx_backward(us i,const vd& x){
  //   vd result(Ns);
  //   assert((i>=0) && (i<=gp));
  //   if (i==0){
  //     result=(data[i+1]()-data[i]())/(x[i+1]-x[i]);
  //     return result;
  //   }
  //   else{
  //     result=(data[i]()-data[i-1]())/(x[i]-x[i-1]);

  //     return result;
  //   }

  // }

  // vd vvar::ddx_central(us i,const vd& x) {
  //   vd result(Ns);
  //   assert((i>=0) && (i<=gp));
  //   if((i>0) && (i<gp-1)) {
  //     result=(data[i+1]()-data[i-1]())/(x[i+1]-x[i-1]); //Central difference
  //     return result;
  //   }
  //   else if(i==0) {
  //     result=(4*data[1]()-3*data[0]()-data[2]())/(2*x[1]);
  //     return result;
  //   }
  //   else {
  //     result=(data[i-2]()-4*data[i-1]()+3*data[i]())/(2*(x[i]-x[i-1]));
  //     return result;
  //   }

  // }
  // us vvar::size() {return gp;}
  // var& vvar::operator[](us i) {return data[i]; } //Extract a variable
  // vvar::~vvar()
  // {
  //   //dtor
  // }





} /* namespace variable */




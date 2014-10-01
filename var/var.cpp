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
    vd newtdata=scalar*var1.tdata();
    return var(var1.gc(),newtdata);
  }

  //***************************************** The var class
  var::var(const Globalconf& gc): var(gc,0.0) {  }
  var::var(const Globalconf& gc1,double initval) :gc__(&gc1),Nf(gc1.Nf()),Ns(gc1.Ns()) {
    TRACE(0,"var::var(const Globalconf& gc, double initval)");
    timedata=vd(gc__->Ns());
    amplitudedata=vd(gc__->Ns());
    settdata(initval);
  }
  var::var(const Globalconf& gc,const vd& timedata):var(gc){ // Create a variable and fill it with time data.
    TRACE(0,"var::var(gc,timedata)");
    this->timedata=timedata;
    dft();
  }
  void var::resetHarmonics(){
    TRACE(10,"var::resetHarmonics()");
    if(gc__->Nf()>1){
      for(us i=3;i<gc__->Ns();i++)
        amplitudedata(i)=0;
    }
    idft();
  }
  var::var(const var& other){
    // Nope, we do not adjust gc__
    this->Nf=other.Nf;
    this->Ns=other.Ns;
    this->timedata=other.timedata;
    this->amplitudedata=other.amplitudedata;
  }
  var& var::operator=(const var& other){
    // THIS WOULD COUPLE TO THE WRONG GLOBALCONF when setRes is used
    // between tubes!!!!!
    if(this->gc__==NULL)
      this->gc__=other.gc__;
    
    this->Nf=other.Nf;
    this->Ns=other.Ns;
    this->timedata=other.timedata;
    this->amplitudedata=other.amplitudedata;
    return *this;
  }
  var var::operator+(const var& other) const{
    var result(*this->gc__);
    result.set(this->operator()()+other());
    return result;
  }
  void var::updateNf(){
    TRACE(5,"var::updateNf(), this Ns:"<< Ns << ", new Ns:" << gc__->Ns() << ". Adress gc__: " << gc__);
    if(this->Ns!=gc__->Ns()){
      assert((gc__->Ns()%2)==1);	// Check if number of samples is not even
      if(this->Ns>gc__->Ns()){
        amplitudedata=amplitudedata.subvec(0,gc__->Ns()-1);
      }
      else{
        vd oldadata=amplitudedata;
        amplitudedata=vd(gc__->Ns(),fillwith::zeros);
        // TRACE(25,"New amplitude data size: "<< amplitudedata.size());
        amplitudedata.subvec(0,this->Ns-1)=oldadata;
      }
      this->Ns=gc__->Ns();		       // Update this number of samples
      timedata=vd(gc__->Ns(),fillwith::zeros); // Reinitialize timedata
      idft();
    }
  }
  var var::operator*(const var& var2) { // Multiply two
    // variables in time domain
    TRACE(0,"var::operator*(const var& var2) const");
    assert(this->Ns==var2.Ns);
    var result(this->gc__);
    vd tdata=this->tdata()%var2.tdata();
    result.settdata(tdata);
    return result;
  }
  var var::operator*(const d& scalar) const {	// Post-multiplication with scalar
    assert(this->gc__!=NULL);
    vd thisadata=this->tdata();
    return var(*(this->gc__),scalar*thisadata);
  }
  
  // Get methods (which require implementation)
  const d& var::operator()(us i) const {//Extract result at specific frequency
    TRACE(-2,"var::operator()("<<i<<"), Ns: "<< Ns);
    assert(i<Ns);
    TRACE(-1,"amplitudedata: "<<amplitudedata);
    return amplitudedata(i);
  }  
  vc var::getcRes() const
  {
    TRACE(-2,"var::getcRes()");
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
    TRACE(-2,"var::tdata(d t)");
    d result=0;
    vc cres=getcRes();
    for(us n=0;n<Nf+1;n++)
      {
        result+=real(cres(n)*exp(I*double(n)*gc__->getomg()*t));
      }
    return result;
  }

  // Set methods
  void var::set(us freqnr,d val) { //Set result for specific frequency zero,real one, -imag one, etc
    TRACE(-2,"var::set("<<freqnr<<","<<val<<")");
    assert(freqnr<Ns);
    amplitudedata[freqnr]=val;
    idft();
    TRACE(-3,"var::set(d val,us freqnr) adata:"<<amplitudedata);
  }
  void var::set(const vc& res)
  {
    TRACE(0,"var::set(const vc& res)");
    assert(res.size()==gc__->Nf()+1);
    amplitudedata(0)=res(0).real();
    for(us i=1;i<Nf+1;i++){
      amplitudedata(2*i-1)=res(i).real();
      amplitudedata(2*i)=res(i).imag();
    }
    idft();
  }
  void var::set(const vd val) {
    TRACE(0,"var::set(const vd& val)");
    assert(val.size()==amplitudedata.size());
    amplitudedata=val;
    idft();
  }
  void var::settdata(double val) {
    TRACE(0,"var::settdata(double val)");
    timedata.fill(val);
    dft();
  }
  void var::settdata(vd& val) {
    TRACE(0,"var::settdata(vd& val)");
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
  dmat var::freqMultiplyMat() const{
    dmat result(gc__->Ns(),gc__->Ns(),fillwith::zeros);
    result(0,0)=amplitudedata(0);
    if(Nf>0){
      for(us j=1;j<gc__->Nf()+1;j++){
        result(2*j-1,2*j-1)= amplitudedata(2*j-1);
        result(2*j-1,2*j  )=-amplitudedata(2*j  ); // Yes only one
        // minus sign
        result(2*j  ,2*j-1)= amplitudedata(2*j);      
        result(2*j  ,2*j  )= amplitudedata(2*j-1);      
      }	// end for loop
    }	// if(Nf>0)
    return result;
  }
  
  // Internal methods for syncing time and amplitude data
  void var::dft() {
    TRACE(0,"var::dft()");
    amplitudedata=gc__->fDFT*timedata;
  }
  void var::idft() { //Internal idft
    timedata=gc__->iDFT*amplitudedata;
  }

  //Get a variable which is the time derivative of the current one
  var var::ddt() const {
    var result(*(this->gc__));
    vd newadata=gc__->DDTfd*amplitudedata;
    result.set(newadata);
    return result;
  }
  // The product

  var var::operator/(const var& var2) const
  {
    vd tdata=this->tdata()/var2.tdata();
    var newvar(*(this->gc__));
    newvar.settdata(tdata);
    return newvar;
  }
  //***************************************** End of the var class



} /* namespace variable */




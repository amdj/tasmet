// var.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////
#include "var.h"
#include "exception.h"
#define Ns (gc_->Ns())
#define Nf (gc_->Nf())
#define fDFT (gc_->fDFT)
#define iDFT (gc_->iDFT)

namespace tasystem {

  using namespace tasystem
;
  //******************************************************************************** Operators
  
  var operator*(const double& scalar,const var& var1){ // Pre-multiplication with scalar
    TRACE(0,"var::operator*(scalar,var)");
    return var(*var1.gc_,scalar*var1.adata_);
  }
  var var::operator/(const var& var2) const {
    TRACE(0,"var::operator/()");
    return var(*gc_,this->tdata()/var2.tdata_,false);
  }
  var var::operator/(const d& val) const {
    TRACE(0,"var::operator/(double)");
    return var(*gc_,this->tdata()/val,false);
  }
  var var::operator-(const var& other) const{
    TRACE(0,"var::operator-(var)");
    return var(*gc_,adata_-other.adata_);
  }
  var var::operator*(d scalar) const {
    TRACE(0,"var::operator*(scalar)");	// Post-multiplication with scalar
    return var(*gc_,scalar*adata_);
  }
  var var::operator+(const var& other) const{
    TRACE(0,"var::operator+()");
    return var(*gc_,adata_+other.adata_);
  }
  var var::operator*(const var& var2) const { // Multiply two
    TRACE(0,"var::operator*(const var& var2) const");
    return var(*gc_,tdata_+var2.tdata_,false);
  }
  //***************************************** The var class
  var::var(const Globalconf& gc): var(gc,0.0) {  }
  var::var(const Globalconf& gc,double initval):
    gc_(&gc)
 {
    TRACE(0,"var::var(const Globalconf& gc, double initval)");
    tdata_=vd(Ns);
    adata_=vd(Ns);
    settdata(initval);
  }
  var::var(const Globalconf& gc,const vd& data,bool adata)
    :var(gc)
  { // Create a tasystem and fill it with time data.
    TRACE(0,"var::var(gc,tdata_)");
    if(data.size()!=adata_.size())
      throw MyError("Wrong size of amplitude vector given. Does the"
                    "vector size correspond to Ns?");
    if(adata){
      this->adata_=data;
      tdata_=iDFT*adata_;
    }
    else{
      this->tdata_=data;
      adata_=fDFT*tdata_;
    }
  }  var::var(const Globalconf& gc,const vc& data)
    :var(gc)
  { // Create a tasystem and fill it with time data.
    if(data.size()!=Nf+1)
      throw MyError("Wrong size of amplitude vector given. Does the"
                    "vector size correspond to Ns?");
    setadata(data);
  }
  void var::setGc(const Globalconf& gc){
    TRACE(10,"var::setGc()");
    this->gc_=&gc;
    updateNf();
  }
  void var::resetHarmonics(){
    TRACE(10,"var::resetHarmonics()");
    if(Nf>1){
      for(us i=3;i<Ns;i++)
        adata_(i)=0;
    }
    tdata_=iDFT*adata_;
  }
  var::var(const var& other){
    TRACE(0,"var::var(const var& other)");
    this->gc_=other.gc_;
    this->tdata_=other.tdata_;
    this->adata_=other.adata_;
  }
  var& var::operator=(const var& other){
    // THIS WOULD COUPLE TO THE WRONG GLOBALCONF when setRes is used
    // between ducts!!!!!
    if(this!=&other){
      if(this->gc_==nullptr)
        this->gc_=other.gc_;
    
      this->tdata_=other.tdata_;
      this->adata_=other.adata_;
      updateNf();
    }
    return *this;
  }
  void var::updateNf(){
    TRACE(0,"var::updateNf()");
    us asize=adata_.size();
    assert(gc_);
    assert(asize==tdata_.size());
    if(asize!=Ns){
      assert((Ns%2)==1);	// Check if number of samples is not even
      if(asize>Ns){
        TRACE(0,"Shrinking vector");
        adata_=adata_.subvec(0,Ns-1);
      }
      else{
        TRACE(0,"Growing vector");
        vd oldadata=adata_;
        adata_=vd(Ns,fillwith::zeros);
        // TRACE(25,"New amplitude data size: "<< adata_.size());
        if(oldadata.size()>0)
          adata_.subvec(0,oldadata.size()-1)=oldadata;
      }
      tdata_=vd(Ns,fillwith::zeros); // Reinitialize tdata_
      tdata_=iDFT*adata_;
    }
  }
  // Get methods (which require implementation)
  d var::operator()(us i) const {//Extract result at specific frequency
    TRACE(-2,"var::operator()("<<i<<"), Ns: "<< Ns);
    if(i>=Ns)
      throw MyError("Invalid frequency number!");
    TRACE(-1,"adata_: "<<adata_);
    return adata_(i);
  }  
  vc var::getcRes() const
  {
    TRACE(-2,"var::getcRes()");
    //	TRACE(0,"adata_:" << adata_);
    vc cadata(Nf+1);
    cadata(0)=adata_(0);
    for(us i=1;i<Nf+1;i++)
      cadata(i)=adata_(2*i-1)+I*adata_(2*i); //The minus is very important
    //	TRACE(0,"Resulting cadata:" << cadata);
    return cadata;
  }
  // Set methods
  void var::setadata(us freqnr,d val) { //Set result for specific frequency zero,real one, -imag one, etc
    TRACE(-2,"var::setadata("<<freqnr<<","<<val<<")");
    assert(freqnr<Ns);
    adata_[freqnr]=val;
    tdata_=iDFT*adata_;
  }
  void var::setadata(const vc& res)
  {
    TRACE(0,"var::setadata(const vc& res)");
    assert(res.size()==Nf+1);
    adata_(0)=res(0).real();
    for(us i=1;i<Nf+1;i++){
      adata_(2*i-1)=res(i).real();
      adata_(2*i)=res(i).imag();
    }
    tdata_=iDFT*adata_;
  }
  d min(us a,us b){return a<=b?a:b;}
  void var::setadata(const vd& val) {
    TRACE(0,"var::setadata(const vd& val)");
    adata_.zeros();
    us minsize=min(val.size(),adata_.size());
    for(us i=0;i<minsize;i++)
      adata_(i)=val(i);
    tdata_=iDFT*adata_;
  }
  void var::settdata(double val) {
    TRACE(0,"var::settdata(double val)");
    for (auto it = tdata_.begin(); it != tdata_.end(); ++it) {
      (*it)=val;
    }
    adata_=fDFT*tdata_;
  }
  void var::settdata(const vd& val) {
    TRACE(0,"var::settdata(vd& val)");
    assert(val.size()==tdata_.size());
    tdata_=val;
    adata_=fDFT*tdata_;
  }

  vd var::timeResponse(us nperiod,us ninst) const {
    TRACE(15,"vd var::timeResponse()");
    vd t=timeResponseTime(nperiod,ninst);
    vc cres=getcRes();
    vd res(t.size(),fillwith::zeros);
    c omg=gc_->getomg();
    for(us i=0;i<Nf+1;i++)
      res+=real(cres(i)*exp(((d) i)*omg*I*t));
    return res;
  }
  vd var::timeResponseTime(us nperiod,us ninst) const{
    d T=1/gc_->getfreq();
    return linspace(0,nperiod*T,ninst);
  }
  dmat var::freqMultiplyMat() const{
    TRACE(0,"var::freqMultiplyMat()");
    dmat result(Ns,Ns,fillwith::zeros);
    result(0,0)=adata_(0);
    if(Nf>0){
      for(us j=1;j<Nf+1;j++){
        result(2*j-1,2*j-1)= adata_(2*j-1);
        result(2*j-1,2*j  )=-adata_(2*j  ); // Yes only one
        // minus sign
        result(2*j  ,2*j-1)= adata_(2*j);      
        result(2*j  ,2*j  )= adata_(2*j-1);      
      }	// end for loop
    }	// if(Nf>0)
    return result;
  }

  //Get a tasystem which is the time derivative of the current one
  var var::ddt() const {
    var result(*(this->gc_));
    vd newadata=gc_->DDTfd*adata_;
    result.setadata(newadata);
    return result;
  }
  // The product

  //***************************************** End of the var class



} /* namespace tasystem */




/*
 * var.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */
#pragma once
#ifndef VAR_H_
#define VAR_H_
#include "globalconf.h"
#include <vtypes.h>
#include <assert.h>


namespace tasystem {
  #ifndef SWIG
  SPOILNAMESPACE
  
  class var;			// Forward declaration
  var operator*(const double&,const var&);
  // var operator+(const d&,const var&);

  #endif
  #ifdef SWIG
  %catches(std::exception,...) var::var(const tasystem::Globalconf&); 
  %catches(std::exception,...) var::var(const tasystem::Globalconf&,double); 
  %catches(std::exception,...) var::var(const tasystem::Globalconf&,const vd& data,bool adata=true); 
  %catches(std::exception,...) var::setadata(const vd& values);
  %catches(std::exception,...) var::settdata(double value);
  %catches(std::exception,...) var::settdata(const vd& values);
  %catches(std::exception,...) var::setadata(us freq,double val);
  %catches(std::exception,...) var::setadata(const vc& values);
  #endif
  class var {
    #ifndef SWIG
    friend var operator*(const double&,const var&);
    #endif
    int dofnr=-1;
    const tasystem::Globalconf* gc_=nullptr;
    vd tdata_,adata_;
  public:
    void setDofNr(us Dofnr){dofnr=Dofnr;}
    int getDofNr() const{return dofnr;}
    var() {}
    var(const var& o);
    var(const tasystem::Globalconf&);	// Initialize with zeros
    var(const tasystem::Globalconf&,double); // Initialize with one time-average value
    var(const tasystem::Globalconf&,const vd& data,bool adata=true); // Initialize with amplitudedata. With tdata_ if adata is set to false
    // Assign with frequency data
    var(const tasystem::Globalconf&,const vc& data);
    

    #ifndef SWIG
    var& operator=(const var&);			  // Copy assignment operator
    #endif
    ~var(){}
    // var operator()(const var&); //Copy constructor
    // Get methods
    const tasystem::Globalconf& gc() const {return *gc_;}
    d operator()(us i) const;				   // Extract amplitude data result at specific frequency    

    // Extract data
    const vd& tdata() const  {return tdata_; } //Get time data
    const vd& adata() const  {return adata_; } //Get time data                                                 //vector

    // Shortcut for adata()
    const vd& operator()() const { return adata_;} //Extract result

    // Obtain a time response vector
    vd timeResponse(us nperiod=2,us ninst=100) const;
    // Obtain the corresponding time vector
    vd timeResponseTime(us nperiod=2,us ninst=100) const;

    #ifndef SWIG
    dmat diagt() const {return diagmat(tdata_);}
    dmat diag() const {return diagmat(adata_);}    
    //Set methods

    void setGc(const tasystem::Globalconf& gc);
    void setGc(const var& o){setGc(*o.gc_);}
    void resetHarmonics();
    void updateNf();

    #endif

    
    void setadata(const vd& values); //Set amplitude data vector to these values
    void settdata(double value); //Set time data to specific value for all time
    void settdata(const vd& values);
    void setadata(us freq,double val); //Set result vector at specific frequency
    void setadata(const vc& values); //Set result vector to these values,

    #ifndef SWIG
    us size() const {return adata_.size();}
    // Specific methods to the result using time domain data
    // Operations ********************

    // Time derivative of this tasystem
    var ddt() const;
    var operator/(const var& var2) const; // Time-domain division operator
    var operator/(const d& var2) const; // Time-domain or frequency domain division operator
    // Multiply two variables in time domain
    var operator*(const var& tasystem) const;
    // Multiply a tasystem with a scalar. This operation is possible
    // for both frequency and time domain data
    var operator*(d scalar) const;

    // add two variables
    var operator+(const var& other) const;
 //Subtract two variables
    var operator-(const var& var2) const;
    // with Note multiplication is defined outside of the class

    // If we need to multiply two numbers in frequency domain, this
    // corresponds to a matrix-vector multiplication (cosines and
    // sines) are mixed up due to the complex numbers. This product
    // can be obtained by getting the matrix-variant of the first
    // tasystem. The following function will give the effective matrix
    dmat freqMultiplyMat() const;
    #endif
    vc getcRes() const;
  };


} /* namespace tasystem */
#endif /* VAR_H_ */


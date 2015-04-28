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


namespace variable {
  #ifndef SWIG
  SPOILNAMESPACE

  class var;			// Forward declaration
  ostream& operator<< (ostream& out,var& v);
  var operator*(const double&,const var&);
  // var operator+(const d&,const var&);

  #endif
  #ifdef SWIG
  %catches(std::exception,...) var::var(const tasystem::Globalconf&); 
  %catches(std::exception,...) var::var(const tasystem::Globalconf&,double); 
  %catches(std::exception,...) var::var(const tasystem::Globalconf&,const vd& data,bool adata=true);
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


    dmat diagt() const {return diagmat(tdata_);}
    dmat diag() const {return diagmat(adata_);}    
    //Set methods

    void setGc(const tasystem::Globalconf& gc);
    void setGc(const var& o){setGc(*o.gc_);}
    void resetHarmonics();
    void updateNf();

    // Extract result at a certain time [s]
    d tdata(d t) const;

    void setadata(const vd& values); //Set amplitude data vector to these values
    void settdata(double value); //Set time data to specific value for all time
    void settdata(const vd& values);
    void setadata(us freq,double val); //Set result vector at specific frequency
    void setadata(const vc& values); //Set result vector to these values,
    us size() const {return adata_.size();}
    // Specific methods to the result using time domain data
    //Show methods
    void showtdata() const; //Print time data to
    void showRes() const;

    // Operations ********************

    // Time derivative of this variable
    var ddt() const;
    var operator/(const var& var2) const; // Time-domain division operator

    // Multiply two variables in time domain
    var operator*(const var& variable) const;
    // Multiply a variable with a scalar. This operation is possible
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
    // variable. The following function will give the effective matrix
    dmat freqMultiplyMat() const;
  private:
    vc getcRes() const;
  };


} /* namespace variable */
#endif /* VAR_H_ */
















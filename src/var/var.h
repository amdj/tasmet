/*
 * var.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */
#pragma once
#ifndef VAR_H_
#define VAR_H_
#include <vtypes.h>
#include <assert.h>
#include "globalconf.h"


namespace variable {
  #ifndef SWIG
  SPOILNAMESPACE

  class var;			// Forward declaration
  ostream& operator<< (ostream& out,var& v);
  var operator*(const double& d1,const var& var2);
  // var operator+(const d&,const var&);

  #endif


  class var {
    int dofnr=-1;
    const tasystem::Globalconf* gc__=NULL;
    vd timedata,amplitudedata;
    us Nf=0,Ns=0;

  public:
    void setDofNr(us Dofnr){dofnr=Dofnr;}
    int getDofNr() const{return dofnr;}
    var() {}
    var(const var& o);
    var(const tasystem::Globalconf&);	// Initialize with zeros
    var(const tasystem::Globalconf *g): var(*g){}
    var(const tasystem::Globalconf&,double); // Initialize with one time-average value
    var(const tasystem::Globalconf&,const vd& timedata); // Initialize with timedata!!!!
    var& operator=(const var&);			  // Copy assignment operator
    ~var(){}
    // var operator()(const var&); //Copy constructor
    // Get methods
    const tasystem::Globalconf& gc() const {return *gc__;}
    const d& operator()(us i) const;				   // Extract amplitude data result at specific frequency    

    const vd& operator()() const { return amplitudedata;} //Extract result
						   //vector
    vc getcRes() const; //Implementation for complex amplitude vector
    const vd& tdata() const  {return timedata; } //Get time data
    const vd& adata() const  {return amplitudedata; } //Get time data                                                 //vector
    d tdata(d t) const; //Extract the estimated value for a given time
                        //t
    dmat diagt() const {return diagmat(timedata);}
    dmat diag() const {return diagmat(amplitudedata);}    
    //Set methods

    void setGc(const tasystem::Globalconf& gc){gc__=&gc;}
    void resetHarmonics();


    void setadata(const vd& values){set(values);}
    #ifndef SWIG
    void set(us freq,double val); //Set result vector at specific frequency
    void set(const vd& values); //Set result vector to these values
    void set(const vc& values); //Set result vector to these values, complex numbers
    #endif
    // Specific methods to the result using time domain data
    void settdata(double value); //Set time data to specific value for all time
    void settdata(const vd& values);
    //Show methods
    void showtdata() const; //Print time data to
    void showRes() const;
    // Operations
    var ddt() const;			  // Time derivative of this variable
    var operator/(const var& var2) const; // Time-domain division operator
    var operator*(const var& variable) const;		   // Multiply two variables in time domain
    var operator*(d scalar) const;   // Multiply a variable with a scalar. This operation is possible for both
				      // frequency and time domain data
    var operator+(const var& other) const;  // add two variables
    var operator-(const var& var2) const; //Subtract two variables
    // with Note multiplication is defined outside of the class

    // If we need to multiply two numbers in frequency domain, this
    // corresponds to a matrix-vector multiplication (cosines and
    // sines) are mixed up due to the complex numbers. This product
    // can be obtained by getting the matrix-variant of the first
    // variable. The following function will give the effective matrix
    dmat freqMultiplyMat() const;
    void updateNf();
    
  protected:
    void dft();
    void idft();

  };


} /* namespace variable */
#endif /* VAR_H_ */
















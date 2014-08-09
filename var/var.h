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
  SPOILNAMESPACE

  using namespace tasystem;
  class var;			// Forward declaration
  ostream& operator<< (ostream& out,var& v);
  var operator*(const double& d1,const var& var2);
  // var operator+(const d&,const var&);

  class var {
  protected:
    vd timedata,amplitudedata;
  public:
    const Globalconf* gc=NULL;
    us Nf=0,Ns=0;
    
    var() {}
    var(const Globalconf&);	// Initialize with zeros
    var(const Globalconf *g): var(*g){}
    var(const Globalconf&,double); // Initialize with one time-average value
    var(const Globalconf&,const vd& timedata); // Initialize with timedata!!!!
    var& operator=(const var&);			  // Copy assignment operator
    // var operator()(const var&); //Copy constructor
    // Get methods

    const d& operator()(us i) const;				   // Extract amplitude data result at specific frequency
    d& operator()(us i);				   // Extract amplitude data result at specific frequency    

    const vd& operator()() const { return amplitudedata;} //Extract result
						   //vector
    var operator*(const var& variable) const;		   // Multiply two variables in time domain
    var operator*(const d& scalar) const;   // Multiply a variable with a scalar. This operation is possible for both
				      // frequency and time domain data
    var operator+(const var& other);  // add two variables

    vd getResfluc() const { return amplitudedata.subvec(1,Ns-1);}
    vc getcRes() const; //Implementation for complex amplitude vector
    const vd& tdata() const  {return timedata; } //Get time data vector
    d tdata(d t) const; //Extract the estimated value for a given time t
    dmat diagt() const {return diagmat(timedata);}
    dmat diag() const {return diagmat(amplitudedata);}    
    //Set methods
    void set(us freq,double val); //Set result vector at specific frequency
    void set(const vd values); //Set result vector to these values
    void set(const vc& values); //Set result vector to these values, complex numbers
    void setResfluc(vd& values); //Set result vector for only unsteady Fourier components
    // Specific methods to the result using time domain data
    void settdata(double value); //Set time data to specific value for all time
    void settdata(vd& values);
    //Show methods
    void showtdata() const; //Print time data to
    void showRes() const;
    // Operations
    var ddt() const;			  // Derivative operator
    var operator/(const var& var2) const; // Time-domain division operator
    var operator-(const var& var2) const; //Not yet implemented
    var operator*(d scalar);			   // post-multiplication
    // with Note multiplication is defined outside of the class

    // If we need to multiply two numbers in frequency domain, this
    // corresponds to a matrix-vector multiplication (cosines and
    // sines) are mixed up due to the complex numbers. This product
    // can be obtained by getting the matrix-variant of the first
    // variable. The following function will give the effective matrix
    dmat freqMultiplyMat() const;

    
  protected:
    void dft();
    void idft();
    void updateNf();
  };


} /* namespace variable */
#endif /* VAR_H_ */
















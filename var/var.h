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

  // void setfreq(d freq);
  // void setNf(us Nf);
  using namespace tasystem;
  class var;			// Forward declaration
  ostream& operator<< (ostream& out,var& v);
  var operator*(const double& d1,const var& var2);


  class var {
  protected:
    vd timedata,amplitudedata;
  public:
    const Globalconf* gc=NULL;
    us Nf=0,Ns=0;
    
    
    var() {}
    var(const Globalconf&);	// Initialize with zeros
    var(const Globalconf&,double); // Initialize with one time-average value
    var(const Globalconf&,const vd& timedata); // Initialize with timedata!!!!
    // var& operator=(const var&);			  // Copy assignment operator
    // var operator()(const var&); //Copy constructor
    // Get methods
    ~var();
    const d& operator()(us i) const;				   // Extract amplitude data result at specific frequency
    d& operator()(us i);				   // Extract amplitude data result at specific frequency    

    vd operator()() const { return amplitudedata;} //Extract result
						   //vector
    var operator*(const var& variable) const;		   // Multiply two variables in time domain
    var operator*(const d& scalar) const;   // Multiply a variable with a scalar. This operation is possible for both
				      // frequency and time domain data
    

    vd getResfluc() const { return amplitudedata.subvec(1,Ns-1);}
    vc getcRes() const; //Implementation for complex amplitude vector
    vd tdata() const  {return timedata; } //Get time data vector
    d tdata(d t) const; //Extract the estimated value for a given time t
    //Set methods
    void set(us freq,double val); //Set result vector at specific frequency
    void set(const vd values); //Set result vector to these values
    void set(const vc values); //Set result vector to these values, complex numbers
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
						   // with
						   // scalar. Note:
						   // pre-multiplication
						   // is defined
						   // outside of the
						   // class
  protected:
    void dft();
    void idft();
    void updateNf();
  };


  // class vvar { //variable container, including some methods to get and set the data simultaneously
  // public:
  //   vvar(string name,us,us,double); //name,gridpoints,number of frequencies, initvals
  //   vvar(string name,us,us); //name,gridpoints,number of frequencies, initvals=0
  //   vvar(); //empty constructor for stack allocation

  //   virtual ~vvar();
  //   vd& getRes(); //Get complete result vector for this variable
  //   vd getRes(us); //Get a result vector for specific frequency
  //   void setRes(vd& res,us freq); //Set result for specific
  //   void setRes(vd& res); //Set complete result vector
  //   vd getTdata(); //Extract a vector with all time data
  //   void setTdata(vd); //Set a complete time domain data vector
  //   void showResult(); //Show result vector

  //   var& operator[](us); //return reference to var at position given
  //   vd ddx_central(us i,const vd& x) ; //Compute ddx of this variable in frequency domain using second order central difference method
  //   vd ddx_forward(us i,const vd& x);
  //   vd ddx_backward(us i,const vd& x);
  //   //vd dotn(us,vd);
  //   us size();
  //   string name();
  //   string Name; //Name string for this variable (density,pressure,etc)
  //   //	dmat ddt;
  //   //	dmat fF,iF; //Forward and inverse FFT matrices.
  //   vd Error,Result; //Error vector for this variable, to r/w store info in
  // protected:
  //   std::vector<var> data;	//The data at the gridpoints, cannot be accessed directly only via []

  //   us gp,Nf,Ns,Dofs; //Degrees of freedom in variable, Number of time samples, Number of frequencies, number of gridpoints

  // private:
  // };



} /* namespace variable */
#endif /* VAR_H_ */
















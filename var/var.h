/*
 * var.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */
#pragma once
#ifndef VAR_H_
#define VAR_H_
#include "../common/vtypes.h"
#include <assert.h>

namespace variable {

  class varoperations;
  class var;

  var operator*(const var& var1,const var& var2);

  class var {
  public:
    var(const varoperations&);
    var(const varoperations&,double);
    var& operator=(const var&);
    var& operator()(const var&); //Copy constructor
    // Get methods
    d operator()(us i) const {//Extract result at specific frequency
      TRACELOG("var::operator(us i)");
      assert(i<Ns);
      return amplitudedata(i); }
    vd operator()() const { return amplitudedata;} //Extract result vector
    vd getResfluc() const { return amplitudedata.subvec(1,Ns-1);}
    vc getcRes() const; //Implementation for complex amplitude vector
    vd tdata() const {return timedata; } //Get time data vector
    d tdata(d t) const; //Extract the estimated value for a given time t
    //Set methods
    void set(double,us); //Set result vector at specific frequency
    void set(const vd& values); //Set result vector to these values
    void set(const vc& values); //Set result vector to these values, complex numbers
    void setResfluc(vd& values); //Set result vector for only unsteady Fourier components
    // Specific methods to the result using time domain data
    void settdata(double value); //Set time data to specific value for all time
    void settdata(vd& values);
    //Show methods
    void showtdata(); //Print time data to
    void showRes();
    // Operations
    var ddt() const;
    var operator/(const var& var2) const;
    var operator-(const var& var2) const; //Not yet implemented
    vd prod(const vd&,const vd&);
    vd dft_copy(const vd&);

    virtual ~var();
    const varoperations* vop; //Pointer to varoperations instance
    const us Nf,Ns;

  protected:
    vd timedata,amplitudedata;
    void dft();
    void idft();



  };

  class varoperations
  {
  public:
    varoperations(us Nf,d freq);
    virtual ~varoperations();
    void setfreq(d newfreq);
    dmat iDFT; //inverse discrete Fourier transform matrix
    dmat fDFT; //forward discrete Fourier transform matrix
    dmat DDTfd;//Derivative in frequency domain
    dmat DDTtd;//Derivative in time domain
    dmat ddt; //Derivative matrix only nonzero frequency components
    dmat iddt; //Inverse of derivative matrix only nonzero frequency components
    d omg;
    us Nf,Ns;
    vd omgvec;
  protected:
    void updateiDFT();
    void updatefDFT();
    void updateiomg();


    d oldomg; //Previous omega
  private:
  };

  ostream& operator<< (ostream& out,var& v);

  class vvar { //variable container, including some methods to get and set the data simultaneously
  public:
    vvar(string name,us,us,double); //name,gridpoints,number of frequencies, initvals
    vvar(string name,us,us); //name,gridpoints,number of frequencies, initvals=0
    vvar(); //empty constructor for stack allocation

    virtual ~vvar();
    vd& getRes(); //Get complete result vector for this variable
    vd getRes(us); //Get a result vector for specific frequency
    void setRes(vd& res,us freq); //Set result for specific
    void setRes(vd& res); //Set complete result vector
    vd getTdata(); //Extract a vector with all time data
    void setTdata(vd); //Set a complete time domain data vector
    void showResult(); //Show result vector

    var& operator[](us); //return reference to var at position given
    vd ddx_central(us i,const vd& x) ; //Compute ddx of this variable in frequency domain using second order central difference method
    vd ddx_forward(us i,const vd& x);
    vd ddx_backward(us i,const vd& x);
    //vd dotn(us,vd);
    us size();
    string name();
    string Name; //Name string for this variable (density,pressure,etc)
    //	dmat ddt;
    //	dmat fF,iF; //Forward and inverse FFT matrices.
    vd Error,Result; //Error vector for this variable, to r/w store info in
  protected:
    std::vector<var> data;	//The data at the gridpoints, cannot be accessed directly only via []

    us gp,Nf,Ns,Dofs; //Degrees of freedom in variable, Number of time samples, Number of frequencies, number of gridpoints

  private:
  };



} /* namespace variable */
#endif /* VAR_H_ */

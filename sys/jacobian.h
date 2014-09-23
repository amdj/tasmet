#pragma once
#ifndef _JACBLOCK_H_
#define _JACBLOCK_H_
#include "var.h"
#include "arma_eigen.h"

namespace tasystem{

  class TripletList;
  
  class JacCol{                // Column block of Jacobian
    bool tobeadded=true;
    us __coldof;                  // First dof of column
    dmat __data;                  // data
  public:
    JacCol(const variable::var&);
    JacCol(const variable::var&,const dmat&);
    bool isToAdd() const {return tobeadded;}
    void setToAdd(bool set){tobeadded=set;} // If this is set to
    // false, this column will not be added to the row.
    JacCol& operator+=(const dmat& add);
    dmat& data(){return __data;}
    const dmat& const_data() const {return __data;}
    us colDof() const{return __coldof;}
  };

  class JacRow{                 // Row in Jacobian matrix
  public:
    us __rowdof;                  // Number of first row
    vector<JacCol> jaccols;     // Column blocks
    JacRow(us rowdofnr,us cols=6);
    JacRow(const variable::var&,us cols=6); // Pick rowdofnr from a variable
    void addCol(const JacCol& jaccol);
    JacRow& operator+=(const JacCol& jaccol){addCol(jaccol); return *this;}
    us rowDof() const {return __rowdof;}
  };

  class Jacobian{

  public:
    vector<JacRow> jacrows;
    void operator+=(const Jacobian&);
    void operator+=(const JacRow&);
    TripletList getTriplets() const;
    
  };

  
} // namespace tasystem

#endif /* _JACBLOCK_H_ */

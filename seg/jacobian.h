#pragma once
#ifndef _JACBLOCK_H_
#define _JACBLOCK_H_
#include "var.h"
#include "arma_eigen.h"

namespace tasystem{

  class JacCol{                // Column block of Jacobian

    us __coldof;                  // First dof of column
    dmat __data;                  // data
  public:
    JacCol(const variable::var&);
    JacCol(const variable::var&,const dmat&);
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
    void addCol(const JacCol& jaccol){jaccols.push_back(jaccol);}
    JacRow& operator+=(const JacCol& jaccol){addCol(jaccol); return *this;}
    us rowDof() const {return __rowdof;}
  };

  class Jacobian{

  public:
    vector<JacRow> jacrows;
    void operator+=(const Jacobian&);
    void operator+=(const JacRow&);
    vtriplet getTriplets() const;
    
  };

  
} // namespace tasystem

#endif /* _JACBLOCK_H_ */

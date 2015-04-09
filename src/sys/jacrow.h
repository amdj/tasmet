// jacrow.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef JACROW_H
#define JACROW_H

#include "jaccol.h"
namespace variable{ class var;}

namespace tasystem{

  class JacRow{                 // Row in Jacobian matrix
    int rowdof_=-1;            // Number of first row, default is
                                // invalid state
  public:
    // Negate all terms
    JacRow operator-() const;

    vector<JacCol> jaccols;     // Column blocks
    JacRow(const JacCol&);
    JacRow(int rowdofnr,us cols=6);
    JacRow(us cols=2):JacRow(-1,cols){}
    void addCol(const JacCol& jaccol);
    JacRow& operator+=(const JacCol& jaccol){addCol(jaccol); return *this;}
    JacRow& operator+=(const JacRow& jacrow);
    JacRow& operator*=(const d& val); // Multiply all terms with constant value

    const int& getRowDof() const {return rowdof_;}
    void setRowDof(us dofnr){rowdof_=dofnr;}
    void show() const;
  };

  
} // namespace tasystem


#endif // JACROW_H
//////////////////////////////////////////////////////////////////////

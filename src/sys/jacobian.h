#pragma once
#ifndef _JACBLOCK_H_
#define _JACBLOCK_H_
#include "arma_eigen.h"

namespace variable{ class var;}
namespace tasystem{ class Globalconf;}

namespace tasystem{
  SPOILNAMESPACE
  class TripletList;
  class JacCol{                // Column block of Jacobian
    bool tobeadded=true;
    us __coldof;                  // First dof of column
    dmat __data;                  // data
  public:
    JacCol(const variable::var&);
    JacCol(us coldof,const tasystem::Globalconf*);
    JacCol(us coldof,const dmat&);    
    JacCol(const variable::var&,const dmat&);
    bool isToAdd() const {return tobeadded;}
    void setToAdd(bool set){tobeadded=set;} // If this is set to
    // false, this column will not be added to the row.
    JacCol& operator+=(const dmat& add);
    dmat& data(){return __data;}
    const dmat& const_data() const {return __data;}
    const us& getColDof() const{return __coldof;}
    void show() const;
  };

  class JacRow{                 // Row in Jacobian matrix
    int __rowdof=-1;            // Number of first row, default is
                                // invalid state
  public:
    vector<JacCol> jaccols;     // Column blocks
    JacRow(int rowdofnr,us cols=6);
    JacRow(us cols):JacRow(-1,cols){}
    JacRow(const variable::var&,us cols=6); // Pick rowdofnr from a variable
    void addCol(const JacCol& jaccol);
    JacRow& operator+=(const JacCol& jaccol){addCol(jaccol); return *this;}
    JacRow& operator+=(const JacRow& jacrow);
    JacRow& operator*=(const d& val); // Multiply all terms with constant value

    const int& getRowDof() const {return __rowdof;}
    void setRowDof(us dofnr){__rowdof=dofnr;}
    void show() const;
  };

  class Jacobian{

  public:
    vector<JacRow> jacrows;
    us maxRow() const;
    us maxCol() const;
    void operator+=(const Jacobian&);
    void operator+=(const JacRow&);
    TripletList getTriplets() const;
    
  };

  
} // namespace tasystem

#endif /* _JACBLOCK_H_ */

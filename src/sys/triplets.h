#pragma once
#ifndef _TRIPLETS_H_
#define _TRIPLETS_H_

#include "vtypes.h"

namespace tasystem{
  SPOILNAMESPACE
  class Triplet{
    int row_,col_;
    d value_;
    bool valid=true;
  public:
    Triplet(const int& row,const int& col,const d& value):row_(row),col_(col),value_(value){}
    const int& col() const {
      // TRACE(25,"row: "<< row()<<", col "<< col_<< ", value:"<< value());
      return col_;
    }
    const int& row() const {return row_;}
    const d& value() const {return value_;}
    void setInvalid(){valid=false;}
    const bool& isValid() const {return valid;}
    void setValue(const d& value){ value_=value;}
    void setRow(const int& r){row_=r;}
    void setCol(const int& c){col_=c;}    
  };



  class TripletList: public vector<Triplet> {
    bool valid=true;
    void setValid_();                                    // Do this if invalid

    // Size of the matrix to build eventually
    us ndofs_;
  public:
    TripletList(us ndofs):ndofs_(ndofs){}
    virtual ~TripletList(){}
    void setValid(){ if(!valid) setValid_();}            // Updates trlist to remove invalid elements

    // Convert to Armadillo Sparse matrix
    operator arma::sp_mat() const;

    void show() const;
    // Make one row zero
    void zeroOutRow(us rownr);
    void multiplyTriplets(const d& multiplicationfactor);

    // Add to capacity
    void reserveExtraDofs(us n);
    
    // Shift position of triplets a certain number of rows and cols.
    void shiftTriplets(int nrows,int ncols);

  };


  





} // namespace tasystem

#endif /* _TRIPLETS_H_ */

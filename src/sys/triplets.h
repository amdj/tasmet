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
    us nrows_,ncols_;
  public:
    TripletList(us nrows,us ncols):nrows_(nrows),ncols_(ncols){}
    virtual ~TripletList(){}
    void setValid(){ if(!valid) setValid_();}            // Updates trlist to remove invalid elements
    // Convert to Armadillo Sparse matrix
    operator arma::sp_mat() const;
    void show() const;
    void zeroOutRow(us rownr);
    void multiplyTriplets(const d& multiplicationfactor);
    void reserveExtraDofs(us n); // Add to capacity
    void shiftTriplets(int nrows,int ncols); // Shift
    // position of triplets a certain number of rows and cols.

  };

  // TripletList getTriplets(const esdmat& mat);
  // TripletList getTripletsBlock(const esdmat& mat,us startrow,us startcol,us nrows,us ncols);

  





} // namespace tasystem

#endif /* _TRIPLETS_H_ */

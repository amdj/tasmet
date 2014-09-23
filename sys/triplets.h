#pragma once
#ifndef _TRIPLETS_H_
#define _TRIPLETS_H_

#include "arma_eigen.h"
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



  class TripletList {
    bool valid=true;
    void setValid_();                                    // Do this if invalid
  public:
    vector<Triplet> trlist;
    TripletList(){}

    virtual ~TripletList(){}
    void setValid(){ if(!valid) setValid_();}            // Updates trlist to remove invalid elements

    void push_back(const Triplet& t){trlist.push_back(t);}
    void reserve(us n) {trlist.reserve(n);}
    int size() const { return trlist.size();}
    Triplet& operator[](int n){return trlist[n];}
    Triplet& at(us n){assert(n<trlist.size()); return trlist[n];}

    vector<Triplet>::iterator begin() { return trlist.begin();}
    vector<Triplet>::iterator end() {return trlist.end();}

    void zeroOutRow(us rownr);
    void multiplyTriplets(const d& multiplicationfactor);
    void reserveExtraDofs(us n); // Add to capacity
    void shiftTriplets(int nrows,int ncols); // Shift
    // position of triplets a certain number of rows and cols.

  };

  TripletList getTriplets(const esdmat& mat);
  TripletList getTripletsBlock(const esdmat& mat,us startrow,us startcol,us nrows,us ncols);

  





} // namespace tasystem

#endif /* _TRIPLETS_H_ */

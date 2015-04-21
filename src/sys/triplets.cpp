#include "triplets.h"

namespace tasystem{
  using arma::sp_mat;

  TripletList::operator sp_mat() const {
    TRACE(15,"TripletList::operator sp_mat()");
    us nvals=size();
    arma::umat locations(2,nvals);
    vd values(nvals);
    for(us i=0;i<nvals;i++) {
      locations(0,i)=(*this)[i].row();
      locations(1,i)=(*this)[i].col();
      values(i)=(*this)[i].value();
    }
    return sp_mat(true,locations,values,ndofs_,ndofs_);
  }
  void TripletList::zeroOutRow(us rownr){
    TRACE(15,"zeroOutRow()");
    valid=false;    
    for(Triplet& tr: *this)
      if(tr.row()==(int) rownr){
        tr.setInvalid();
      }
  }

  void TripletList::setValid_(){ // Set valid if invalid
    TRACE(15,"TripletList::setValid()");
    vector<Triplet> validtrlist;
    reserve(size());
    for(auto tr=begin();tr!=end();++tr)
      if(tr->isValid()){
        validtrlist.push_back(*tr);
      }
    clear();
    vector<Triplet>::operator=(validtrlist);

    valid=true;
  }
  void TripletList::show()  const {
    for(const auto& t: *this){
      cout << "Row: " << t.row() << " , column: " << t.col() << " , value: " << t.value() << "\n";
    }
  }
  void TripletList::multiplyTriplets(const d& factor){
    TRACE(15,"multiplyTriplets()");
    for(auto tr=begin();tr!=end();tr++)
      tr->setValue(factor*tr->value());
  }

  void TripletList::reserveExtraDofs(us n){
    TRACE(15,"reserveExtraDofs()");
    us cursize=size();
    reserve(cursize+n);
  }
  void TripletList::shiftTriplets(int nrows,int ncols){
    // shift the position of the values in a matrix. nrows and ncols
    // can be negative numbers.
    TRACE(15,"shiftTriplets()");
    for(auto tr=begin();tr!=end();tr++){
      tr->setCol(tr->col()+ncols);
      tr->setRow(tr->row()+nrows);
    }
  }


  // TripletList getTripletsBlock(const esdmat& mat,us startrow,us startcol,us nrows,us ncols){
  //   assert(startrow+nrows <= (us) mat.rows());
  //   assert(startcol+ncols <= (us) mat.cols());
  //   us Mj,Mi,i,j,currOuterIndex,nextOuterIndex;
  //   TripletList tripletList;
  //   tripletList.reserve(mat.nonZeros());

  //   for(j=0; j<ncols; j++){
  //     Mj=j+startcol;
  //     currOuterIndex = mat.outerIndexPtr()[Mj];
  //     nextOuterIndex = mat.outerIndexPtr()[Mj+1];

  //     for(us a = currOuterIndex; a<nextOuterIndex; a++){
  //       Mi=mat.innerIndexPtr()[a];

  //       if(Mi < startrow) continue;
  //       if(Mi >= startrow + nrows) break;

  //       i=Mi-startrow;    
  //       tripletList.push_back(Triplet(i,j,mat.valuePtr()[a]));
  //     }
  //   }
  //   return tripletList;
  // }


}            // namespace tasystem

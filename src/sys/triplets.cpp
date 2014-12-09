#include "triplets.h"

namespace tasystem{

  void TripletList::zeroOutRow(us rownr){
    TRACE(15,"zeroOutRow()");
    valid=false;    
    for(auto tr=trlist.begin();tr!=trlist.end();tr++)
      if(tr->row()==(int) rownr){
        tr->setInvalid();
      }
  }

  void TripletList::setValid_(){ // Set valid if invalid
    TRACE(15,"TripletList::setValid()");
    vector<Triplet> validtrlist;
    validtrlist.reserve(trlist.size());
    for(auto tr=trlist.begin();tr!=trlist.end();++tr)
      if(tr->isValid()){
        validtrlist.push_back(*tr);
      }
    trlist.clear();
    trlist=validtrlist;

    valid=true;
  }
  void TripletList::show()  const {
    for(auto t=trlist.begin();t!=trlist.end();t++){
      cout << "Row: " << t->row() << " , column: " << t->col() << " , value: " << t->value() << "\n";
    }
  }
  void TripletList::multiplyTriplets(const d& factor){
    TRACE(15,"multiplyTriplets()");
    for(auto tr=trlist.begin();tr!=trlist.end();tr++)
      tr->setValue(factor*tr->value());
  }

  void TripletList::reserveExtraDofs(us n){
    TRACE(15,"reserveExtraDofs()");
    us cursize=trlist.size();
    trlist.reserve(cursize+n);
  }
  void TripletList::shiftTriplets(int nrows,int ncols){
    // shift the position of the values in a matrix. nrows and ncols
    // can be negative numbers.
    TRACE(15,"shiftTriplets()");
    for(auto tr=trlist.begin();tr!=trlist.end();tr++){
      tr->setCol(tr->col()+ncols);
      tr->setRow(tr->row()+nrows);
    }
  }

  TripletList getTriplets(const esdmat & mat){
    //only for ColMajor Sparse Matrix
    TripletList tripletlist;
    tripletlist.reserve(mat.nonZeros());
    for (int k=0; k<mat.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          it.value();
          it.row(); // row index
          it.col(); // col index (here it is equal to k)
          it.index(); // inner index, here it is equal to it.row()
          tripletlist.push_back(Triplet(it.row(),it.col(),it.value()));
        }
    }
    
    return tripletlist;
  } // getTriplets

  TripletList getTripletsBlock(const esdmat& mat,us startrow,us startcol,us nrows,us ncols){
    assert(startrow+nrows <= (us) mat.rows());
    assert(startcol+ncols <= (us) mat.cols());
    us Mj,Mi,i,j,currOuterIndex,nextOuterIndex;
    TripletList tripletList;
    tripletList.reserve(mat.nonZeros());

    for(j=0; j<ncols; j++){
      Mj=j+startcol;
      currOuterIndex = mat.outerIndexPtr()[Mj];
      nextOuterIndex = mat.outerIndexPtr()[Mj+1];

      for(us a = currOuterIndex; a<nextOuterIndex; a++){
        Mi=mat.innerIndexPtr()[a];

        if(Mi < startrow) continue;
        if(Mi >= startrow + nrows) break;

        i=Mi-startrow;    
        tripletList.push_back(Triplet(i,j,mat.valuePtr()[a]));
      }
    }
    return tripletList;
  }




}            // namespace tasystem

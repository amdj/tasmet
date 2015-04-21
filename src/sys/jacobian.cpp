#include "jacobian.h"
#include "triplets.h"

namespace tasystem{

  void Jacobian::operator+=(const Jacobian& other){
    TRACE(2,"Jacobian::append()");
    jacrows.insert(jacrows.end(),other.jacrows.begin(),other.jacrows.end());
  }
  void Jacobian::operator+=(const JacRow& other){
    TRACE(2,"Jacobian::append()");
    jacrows.push_back(other);
  }
  Jacobian::operator TripletList() const{
    TRACE(18,"Jacobian::operator Tripletlist()");
    int insertrow,insertcol;
    TripletList res(ndofs_);
    const dmat& typicaldatacel=jacrows.at(0).jaccols.at(0).const_data();
    us size=typicaldatacel.n_rows;

    res.reserve(4*jacrows.size()*pow(size,2)); // Should
    // be approximately enough
    us i,j;

    for(const JacRow& row: jacrows) {
      insertrow=row.getRowDof();
      for(const JacCol& col: row.jaccols){
        insertcol=col.getColDof();
        if(insertcol>=0){
          const dmat& data=col.const_data();
          for(i=0;i<size;i++){
            for(j=0;j<size;j++){
              if(data(i,j)!=0)
                res.push_back(Triplet(i+insertrow,j+insertcol,data(i,j)));
            }
          }
        } // insertcol>0
      }   // for loop over cols
    }     // for loop over rows
    // TRACE(10,"SFSG");
    return res;
  }
  
} // namespace tasystem

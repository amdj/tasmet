#include "jacobian.h"

namespace tasystem{

  JacCol::JacCol(const variable::var& thevar):
    __coldof(thevar.getDofNr()),
    __data(thevar.gc->Ns,thevar.gc->Ns,fillwith::zeros)
  {  }
  JacCol::JacCol(const variable::var& thevar,const dmat& data):
    __coldof(thevar.getDofNr()),
    __data(data)
  {  }
  
  JacCol& JacCol::operator+=(const dmat& data){
    __data+=data;
    return *this;
  }

  JacRow::JacRow(const variable::var& vardof,us nrofcols):JacRow(vardof.getDofNr(),nrofcols){}
  JacRow::JacRow(us rowdofnr,us nrofcols): __rowdof(rowdofnr){
    jaccols.reserve(nrofcols);
  }  
  void JacRow::addCol(const JacCol& jaccol){
    if(jaccol.isToAdd())
      jaccols.push_back(jaccol);
  }  

  void Jacobian::operator+=(const Jacobian& other){
    TRACE(2,"Jacobian::append()");
    jacrows.insert(jacrows.end(),other.jacrows.begin(),other.jacrows.end());
  }
  void Jacobian::operator+=(const JacRow& other){
    TRACE(2,"Jacobian::append()");
    jacrows.push_back(other);
  }
  vtriplet Jacobian::getTriplets() const{
    TRACE(18,"Jacobian::getTriplets()");
    us insertrow,insertcol;
    vtriplet res;
    const dmat& typicaldatacel=jacrows.at(0).jaccols.at(0).const_data();
    us size=typicaldatacel.n_rows;

    res.reserve(4*jacrows.size()*pow(size,2)); // Should
    // be approximately enough
    us i,j;
    for(auto row=jacrows.begin();row!=jacrows.end();row++){
      insertrow=row->rowDof();
      for(auto col=row->jaccols.begin();col<row->jaccols.end();col++){
        insertcol=col->colDof();
        const dmat& data=col->const_data();
        for(us i=0;i<size;i++)
          for(us j=0;j<size;j++)
            if(data(i,j)!=0)
              res.push_back(triplet(i+insertrow,j+insertcol,data(i,j)));

      }
    }

    return res;
  }
  
} // namespace tasystem

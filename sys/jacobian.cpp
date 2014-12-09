#include "jacobian.h"
#include "triplets.h"
#include "var.h"
#include "globalconf.h"

namespace tasystem{

  JacCol::JacCol(const variable::var& thevar):
    __coldof(thevar.getDofNr()),
    __data(thevar.gc().Ns(),thevar.gc().Ns(),fillwith::zeros)
  {  }
  JacCol::JacCol(us coldof,const tasystem::Globalconf* gc):
    __coldof(coldof),
    __data(gc->Ns(),gc->Ns(),fillwith::zeros)
  {  }
  JacCol::JacCol(us coldof,const dmat& data):
    __coldof(coldof),
    __data(data)
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
  JacRow::JacRow(int rowdofnr,us nrofcols): __rowdof(rowdofnr){
    jaccols.reserve(nrofcols);
  }  
  void JacRow::addCol(const JacCol& jaccol){
    if(jaccol.isToAdd())
      jaccols.push_back(jaccol);
  }  
  JacRow& JacRow::operator*=(const d& val){
    TRACE(2,"Jacobian::operator*=()");
    for(auto col=jaccols.begin();col!=jaccols.end();col++)
      col->data()*=val;
    return *this;
  }
  JacRow& JacRow::operator+=(const JacRow& other){
    TRACE(2,"Jacobian::operator*=()");
    this->jaccols.reserve(this->jaccols.capacity()+other.jaccols.size()-this->jaccols.size());
    for(auto col=other.jaccols.begin();col!=other.jaccols.end();col++)
      this->operator+=(*col);
    return *this;
  }

  
  
  void Jacobian::operator+=(const Jacobian& other){
    TRACE(2,"Jacobian::append()");
    jacrows.insert(jacrows.end(),other.jacrows.begin(),other.jacrows.end());
  }
  void Jacobian::operator+=(const JacRow& other){
    TRACE(2,"Jacobian::append()");
    jacrows.push_back(other);
  }
  TripletList Jacobian::getTriplets() const{
    TRACE(18,"Jacobian::getTriplets()");
    int insertrow,insertcol;
    TripletList res;
    const dmat& typicaldatacel=jacrows.at(0).jaccols.at(0).const_data();
    us size=typicaldatacel.n_rows;

    res.reserve(4*jacrows.size()*pow(size,2)); // Should
    // be approximately enough
    us i,j;
    // WARN("Dangerous setting");

    for(auto row=jacrows.begin();row!=jacrows.end();row++){
      insertrow=row->getRowDof();
      for(auto col=row->jaccols.begin();col<row->jaccols.end();col++){
        insertcol=col->getColDof();
        if(insertcol>=0){
          const dmat& data=col->const_data();
          for(us i=0;i<size;i++){
            for(us j=0;j<size;j++){
              if(data(i,j)!=0)
            // TRACE(20,"abs data:" << std::abs(data(i,j)));
            // if(std::abs(data(i,j))>1e-15)
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

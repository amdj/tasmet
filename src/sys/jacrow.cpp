#include "jacobian.h"
#include "var.h"
#include "globalconf.h"

namespace tasystem{

  JacRow::JacRow(int rowdofnr,us nrofcols): rowdof_(rowdofnr){
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
  JacRow JacRow::operator-() const{
    TRACE(15,"JacRow::operator-()");
    JacRow result(*this);
    for (auto jaccol = result.jaccols.begin(); jaccol != result.jaccols.end(); ++jaccol) {
      (*jaccol)*=-1;
    }
    return result;
  }
  
} // namespace tasystem

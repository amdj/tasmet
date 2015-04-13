#include "jacobian.h"
#include "var.h"
#include "globalconf.h"

namespace tasystem{

  JacRow::JacRow(const JacCol& j):
    JacRow(-1,1)
  {
    (*this)+=j;
  }
  // JacRow& JacRow::operator+=(JacCol&& j){
  //   TRACE(45,"JacRow::operator+=(JacCol&& j)");
  //   jaccols.emplace_back(std::move(j));
  // }
  JacRow& JacRow::operator+=(const JacCol& j){
    TRACE(10,"JacRow::operator+=(const JacCol& j)");
    if(j.isToAdd())
      jaccols.emplace_back(j);
    return *this;
  }
    
  // JacRow& JacRow::addCol(const JacCol& jaccol){
  //   if(jaccol.isToAdd())
  //     jaccols.push_back(jaccol);
  //   return *this;
  // }  
  JacRow& JacRow::operator*=(const d& val){
    TRACE(2,"Jacobian::operator*=()");
    for(JacCol& col: jaccols)
      col.data()*=val;
    return *this;
  }
  JacRow& JacRow::operator+=(const JacRow& other){
    TRACE(2,"Jacobian::operator*=()");
    this->jaccols.reserve(this->jaccols.capacity()+other.jaccols.size()-this->jaccols.size());
    for(const JacCol& col :other.jaccols)
      (*this)+=col;
    return *this;
  }
  JacRow JacRow::operator-() const{
    TRACE(15,"JacRow::operator-()");
    JacRow result(*this);
    for (JacCol& jaccol : result.jaccols) 
      jaccol*=-1;

    return result;
  }
  void JacRow::prePostMultiply(const dmat& pre,const dmat& post) {
    for(JacCol& jaccol: jaccols)
      jaccol.prePostMultiply(pre,post);
  }
} // namespace tasystem

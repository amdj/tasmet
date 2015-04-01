// jaccol.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "jaccol.h"
#include "var.h"

namespace tasystem{

  JacCol::JacCol(const variable::var& thevar):
    coldof_(thevar.getDofNr()),
    data_(thevar.gc().Ns(),thevar.gc().Ns(),fillwith::zeros)
  {  }
  JacCol::JacCol(us coldof,const tasystem::Globalconf* gc):
    coldof_(coldof),
    data_(gc->Ns(),gc->Ns(),fillwith::zeros)
  {  }
  JacCol::JacCol(us coldof,const dmat& data):
    coldof_(coldof),
    data_(data)
  {  }

  JacCol::JacCol(const variable::var& thevar,const dmat& data):
    coldof_(thevar.getDofNr()),
    data_(data)
  {  }
  
  JacCol& JacCol::operator+=(const dmat& data){
    data_+=data;
    return *this;
  }
  JacCol& JacCol::operator*=(const d& val){
    data_*=val;
    return *this;
  }

} // namespace tasystem

//////////////////////////////////////////////////////////////////////

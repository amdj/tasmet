#pragma once
#ifndef _RIGHTCELL_H_
#define _RIGHTCELL_H_
#include "bccell.h"

namespace tube{

  class RightCell:public BcCell{
    tasystem::var mr_;
  public:
    RightCell(us i,const Tube& t);
    virtual ~RightCell(){}
    virtual void init(const Cell* left,const Cell* right);
    virtual Pos getPos() const {return Pos::right;}
    virtual void show(us detailnr=1) const;

    const d& Sfbc() const {return Sfr;}
    const tasystem::var& mbc() const{return mr_;}
    const tasystem::var& mr() const { return mr_;}
    const tasystem::var& TR() const { return Tbc_;}    


  };  

} // namespace tube

#endif /* _RIGHTCELL_H_ */

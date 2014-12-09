#pragma once
#ifndef _TMTUBES_H_
#define _TMTUBES_H_
#include <wx/wx.h>
#include "globalconf.h"

class Tmtubes: public wxApp
{
public:
  tasystem::Globalconf gc;
  
  virtual bool OnInit();

};

#endif /* _TMTUBES_H_ */

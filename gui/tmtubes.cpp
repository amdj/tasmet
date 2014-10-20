#include "tmtubes.h"
#include "gui.h"

bool Tmtubes::OnInit(){
  gc=tasystem::Globalconf::airSTP(0,100);
  inittrace(15);
  Gui *g=new Gui(this);
  g->Show(true);
    
  return true;
}



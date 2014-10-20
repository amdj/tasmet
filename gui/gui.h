#pragma once
#ifndef _MYTMTUBESGUI_H_
#define _MYTMTUBESGUI_H_
#include <wx/wx.h>
#include "guibase.h"
#include "globalconf.h"

class Tmtubes;

class Gui : public GuiBase
{
  Tmtubes* tmtubes;
public:
  Gui(Tmtubes* tmtubes):GuiBase(NULL),tmtubes(tmtubes){}
protected:
  virtual void editGc( wxCommandEvent& event );

};


class MyGlobalconfDialog : public GlobalconfDialog {
  tasystem::Globalconf* orig=NULL;
  tasystem::Globalconf localcopy;
private:
  void updategc();
  void updatefields();
public:
  MyGlobalconfDialog(wxWindow* parent,Tmtubes* tmtubes);
  // Virtual event handlers, overide them in your derived class
  virtual void updategc( wxKeyEvent& event ) { updategc(); updatefields(); }
  virtual void updategc( wxFocusEvent& event ) {updategc();  updatefields(); }
  virtual void updategc( wxMouseEvent& event ) {updategc();  updatefields(); }
  virtual void OK( wxCommandEvent& event );
  virtual void cancel( wxCommandEvent& event );

};
#endif /* _MYTMTUBESGUI_H_ */


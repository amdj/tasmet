#include "gui.h"
#include "tmtubes.h"
typedef double d;
typedef unsigned us;
inline wxString toString(double d){ return wxString::Format(wxT("%f"),d);}
inline wxString toString(us f){ return wxString::Format(wxT("%d"),f);}



void Gui::editGc( wxCommandEvent& event ){

  MyGlobalconfDialog* gcd=new MyGlobalconfDialog(this,tmtubes);
  gcd->Show(true);
  
}

void errorDialog(const wxString& var){
  wxString msg=wxString::Format("Error parsing variable: %s",var);
  wxMessageDialog *dial = new wxMessageDialog(NULL, 
     msg, wxT("Error"), wxOK | wxICON_ERROR);
   dial->ShowModal();
}

MyGlobalconfDialog::MyGlobalconfDialog(wxWindow* parent,Tmtubes* tmtubes):
  GlobalconfDialog(parent),
  orig(&tmtubes->gc)
{
  localcopy=*orig;
  updatefields();
}

void MyGlobalconfDialog::updatefields(){
  if(localcopy.getGas().compare("air")==0)
    Gas->SetSelection(0);
  if(localcopy.getGas().compare("helium")==0)
    Gas->SetSelection(1);
    
  p0->SetValue(toString(localcopy.p0));
  T0->SetValue(toString(localcopy.T0));
  M->SetValue(toString(localcopy.getMass()));
  freq->SetValue(toString(localcopy.getfreq()));
  Nf->SetValue(toString(localcopy.Nf()));

  rho0->SetValue(toString(localcopy.rho0()));
  c0->SetValue(toString(localcopy.c0()));

}

void MyGlobalconfDialog::updategc(){
  switch(Gas->GetSelection()){
  case 0:
    std::cout << "air selected\n";
    localcopy.setGas("air");
    break;
  case 1:
    std::cout << "helium selected\n";
    localcopy.setGas("helium");
    break;
  }
  
  if(!p0->GetValue().ToDouble(&localcopy.p0)){
    errorDialog(wxT("p0"));
    return;
  }
  if(!T0->GetValue().ToDouble(&localcopy.T0)){
    errorDialog(wxT("T0"));
    return;
  }
  d M1;
  if(!M->GetValue().ToDouble(&M1)){
    errorDialog(wxT("M"));
    return;
  }
  localcopy.setMass(M1);
  d freq1;
  if(!freq->GetValue().ToDouble(&freq1)){
    errorDialog(wxT("f"));
    return;
  }
  localcopy.setfreq(freq1);  
  long Nf1;
  if(!Nf->GetValue().ToLong(&Nf1)){
    errorDialog(wxT("Nf"));
    return;
  }
  localcopy.setNf((us) Nf1);


}
void MyGlobalconfDialog::OK( wxCommandEvent& event ){

  updategc();

  *orig=localcopy;
  Close(true);
}

void MyGlobalconfDialog::cancel( wxCommandEvent& event ){

  Close(true);
}

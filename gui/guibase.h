///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Sep  8 2010)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#ifndef __guibase__
#define __guibase__

#include <wx/string.h>
#include <wx/button.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/sizer.h>
#include <wx/frame.h>
#include <wx/stattext.h>
#include <wx/choice.h>
#include <wx/textctrl.h>
#include <wx/statbox.h>
#include <wx/dialog.h>

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class GuiBase
///////////////////////////////////////////////////////////////////////////////
class GuiBase : public wxFrame 
{
	private:
	
	protected:
		wxButton* edit_globalconf_but;
		
		// Virtual event handlers, overide them in your derived class
		virtual void editGc( wxCommandEvent& event ) { event.Skip(); }
		
	
	public:
		
		GuiBase( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("TMtubes nonlinear"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 500,300 ), long style = wxDEFAULT_FRAME_STYLE|wxTAB_TRAVERSAL );
		~GuiBase();
	
};

///////////////////////////////////////////////////////////////////////////////
/// Class GlobalconfDialog
///////////////////////////////////////////////////////////////////////////////
class GlobalconfDialog : public wxDialog 
{
	private:
	
	protected:
		wxStaticText* m_staticText3;
		wxChoice* Gas;
		wxStaticText* m_staticText5;
		wxStaticText* m_staticText4;
		wxTextCtrl* p0;
		wxStaticText* m_staticText6;
		wxStaticText* m_staticText7;
		wxTextCtrl* T0;
		wxStaticText* m_staticText9;
		wxStaticText* m_staticText10;
		wxTextCtrl* M;
		wxStaticText* m_staticText11;
		wxStaticText* m_staticText12;
		wxTextCtrl* freq;
		wxStaticText* m_staticText13;
		wxStaticText* m_staticText14;
		wxTextCtrl* Nf;
		wxStaticText* m_staticText15;
		wxStaticText* m_staticText17;
		wxChoice* driven;
		wxStaticText* m_staticText18;
		wxTextCtrl* rho0;
		wxStaticText* m_staticText22;
		wxStaticText* m_staticText20;
		wxTextCtrl* c0;
		wxStaticText* m_staticText23;
		wxButton* okbutton;
		wxButton* m_button6;
		
		// Virtual event handlers, overide them in your derived class
		virtual void updategc( wxKeyEvent& event ) { event.Skip(); }
		virtual void updategc( wxFocusEvent& event ) { event.Skip(); }
		virtual void updategc( wxMouseEvent& event ) { event.Skip(); }
		virtual void OK( wxCommandEvent& event ) { event.Skip(); }
		virtual void cancel( wxCommandEvent& event ) { event.Skip(); }
		
	
	public:
		
		GlobalconfDialog( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Global configuration"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 646,530 ), long style = wxDEFAULT_DIALOG_STYLE );
		~GlobalconfDialog();
	
};

#endif //__guibase__

///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Sep  8 2010)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#include "guibase.h"

///////////////////////////////////////////////////////////////////////////

GuiBase::GuiBase( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxFrame( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );
	
	wxBoxSizer* bSizer8;
	bSizer8 = new wxBoxSizer( wxVERTICAL );
	
	edit_globalconf_but = new wxButton( this, wxID_ANY, wxT("Edit global configuration..."), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer8->Add( edit_globalconf_but, 0, wxALL, 5 );
	
	this->SetSizer( bSizer8 );
	this->Layout();
	
	this->Centre( wxBOTH );
	
	// Connect Events
	edit_globalconf_but->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( GuiBase::editGc ), NULL, this );
}

GuiBase::~GuiBase()
{
	// Disconnect Events
	edit_globalconf_but->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( GuiBase::editGc ), NULL, this );
	
}

GlobalconfDialog::GlobalconfDialog( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxDialog( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );
	this->SetBackgroundColour( wxSystemSettings::GetColour( wxSYS_COLOUR_GRAYTEXT ) );
	
	wxBoxSizer* bSizer2;
	bSizer2 = new wxBoxSizer( wxVERTICAL );
	
	wxStaticBoxSizer* sbSizer7;
	sbSizer7 = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Input") ), wxVERTICAL );
	
	wxGridSizer* input;
	input = new wxGridSizer( 8, 3, 0, 0 );
	
	m_staticText3 = new wxStaticText( this, wxID_ANY, wxT("Gas type (g)"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText3->Wrap( -1 );
	input->Add( m_staticText3, 0, wxALL, 5 );
	
	wxString GasChoices[] = { wxT("Air"), wxT("Helium") };
	int GasNChoices = sizeof( GasChoices ) / sizeof( wxString );
	Gas = new wxChoice( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, GasNChoices, GasChoices, 0 );
	Gas->SetSelection( 0 );
	input->Add( Gas, 0, wxALL, 5 );
	
	m_staticText5 = new wxStaticText( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText5->Wrap( -1 );
	input->Add( m_staticText5, 0, wxALL, 5 );
	
	m_staticText4 = new wxStaticText( this, wxID_ANY, wxT("Reference pressure (p0)"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText4->Wrap( -1 );
	input->Add( m_staticText4, 0, wxALL, 5 );
	
	p0 = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 200,-1 ), wxTE_RIGHT );
	input->Add( p0, 0, wxALL, 5 );
	
	m_staticText6 = new wxStaticText( this, wxID_ANY, wxT("Pa"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText6->Wrap( -1 );
	input->Add( m_staticText6, 0, wxALL, 5 );
	
	m_staticText7 = new wxStaticText( this, wxID_ANY, wxT("Reference temperature (T0)"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText7->Wrap( -1 );
	input->Add( m_staticText7, 0, wxALL, 5 );
	
	T0 = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 200,-1 ), wxTE_RIGHT );
	input->Add( T0, 0, wxALL, 5 );
	
	m_staticText9 = new wxStaticText( this, wxID_ANY, wxT("K"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText9->Wrap( -1 );
	input->Add( m_staticText9, 0, wxALL, 5 );
	
	m_staticText10 = new wxStaticText( this, wxID_ANY, wxT("Mass in system (m)"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText10->Wrap( -1 );
	input->Add( m_staticText10, 0, wxALL, 5 );
	
	M = new wxTextCtrl( this, wxID_ANY, wxT("0.0"), wxDefaultPosition, wxSize( 200,-1 ), wxTE_RIGHT );
	input->Add( M, 0, wxALL, 5 );
	
	m_staticText11 = new wxStaticText( this, wxID_ANY, wxT("kg"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText11->Wrap( -1 );
	input->Add( m_staticText11, 0, wxALL, 5 );
	
	m_staticText12 = new wxStaticText( this, wxID_ANY, wxT("Frequency (f)"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText12->Wrap( -1 );
	input->Add( m_staticText12, 0, wxALL, 5 );
	
	freq = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 200,-1 ), wxTE_RIGHT );
	input->Add( freq, 0, wxALL, 5 );
	
	m_staticText13 = new wxStaticText( this, wxID_ANY, wxT("Hz"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText13->Wrap( -1 );
	input->Add( m_staticText13, 0, wxALL, 5 );
	
	m_staticText14 = new wxStaticText( this, wxID_ANY, wxT("Number of harmonics (Nf)"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText14->Wrap( -1 );
	input->Add( m_staticText14, 0, wxALL, 5 );
	
	Nf = new wxTextCtrl( this, wxID_ANY, wxT("0"), wxDefaultPosition, wxDefaultSize, wxTE_RIGHT );
	input->Add( Nf, 0, wxALL, 5 );
	
	m_staticText15 = new wxStaticText( this, wxID_ANY, wxT("-"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText15->Wrap( -1 );
	input->Add( m_staticText15, 0, wxALL, 5 );
	
	m_staticText17 = new wxStaticText( this, wxID_ANY, wxT("Driven system (d)"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText17->Wrap( -1 );
	input->Add( m_staticText17, 0, wxALL, 5 );
	
	wxString drivenChoices[] = { wxT("Yes"), wxT("No") };
	int drivenNChoices = sizeof( drivenChoices ) / sizeof( wxString );
	driven = new wxChoice( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, drivenNChoices, drivenChoices, 0 );
	driven->SetSelection( 0 );
	input->Add( driven, 0, wxALL, 5 );
	
	sbSizer7->Add( input, 1, 0, 5 );
	
	bSizer2->Add( sbSizer7, 1, wxEXPAND, 5 );
	
	wxStaticBoxSizer* computed;
	computed = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Computed values") ), wxHORIZONTAL );
	
	wxGridSizer* gSizer2;
	gSizer2 = new wxGridSizer( 2, 3, 0, 0 );
	
	m_staticText18 = new wxStaticText( this, wxID_ANY, wxT("Reference density (rho0)"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText18->Wrap( -1 );
	gSizer2->Add( m_staticText18, 0, wxALL, 5 );
	
	rho0 = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 200,-1 ), 0 );
	rho0->Enable( false );
	
	gSizer2->Add( rho0, 0, wxALL, 5 );
	
	m_staticText22 = new wxStaticText( this, wxID_ANY, wxT("kg/m^3"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText22->Wrap( -1 );
	gSizer2->Add( m_staticText22, 0, wxALL, 5 );
	
	m_staticText20 = new wxStaticText( this, wxID_ANY, wxT("Referene speed of sound (c0)"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText20->Wrap( -1 );
	gSizer2->Add( m_staticText20, 0, wxALL, 5 );
	
	c0 = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 200,-1 ), 0 );
	c0->Enable( false );
	
	gSizer2->Add( c0, 0, wxALL, 5 );
	
	m_staticText23 = new wxStaticText( this, wxID_ANY, wxT("m/s"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText23->Wrap( -1 );
	gSizer2->Add( m_staticText23, 0, wxALL, 5 );
	
	computed->Add( gSizer2, 1, 0, 5 );
	
	bSizer2->Add( computed, 0, wxEXPAND, 5 );
	
	wxGridSizer* gSizer3;
	gSizer3 = new wxGridSizer( 1, 2, 0, 0 );
	
	okbutton = new wxButton( this, wxID_ANY, wxT("OK"), wxDefaultPosition, wxDefaultSize, 0 );
	gSizer3->Add( okbutton, 0, wxALL, 5 );
	
	m_button6 = new wxButton( this, wxID_ANY, wxT("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
	gSizer3->Add( m_button6, 0, wxALL, 5 );
	
	bSizer2->Add( gSizer3, 1, wxALIGN_CENTER, 5 );
	
	this->SetSizer( bSizer2 );
	this->Layout();
	
	this->Centre( wxBOTH );
	
	// Connect Events
	Gas->Connect( wxEVT_CHAR, wxKeyEventHandler( GlobalconfDialog::updategc ), NULL, this );
	Gas->Connect( wxEVT_KEY_DOWN, wxKeyEventHandler( GlobalconfDialog::updategc ), NULL, this );
	p0->Connect( wxEVT_KILL_FOCUS, wxFocusEventHandler( GlobalconfDialog::updategc ), NULL, this );
	p0->Connect( wxEVT_LEAVE_WINDOW, wxMouseEventHandler( GlobalconfDialog::updategc ), NULL, this );
	T0->Connect( wxEVT_KILL_FOCUS, wxFocusEventHandler( GlobalconfDialog::updategc ), NULL, this );
	freq->Connect( wxEVT_KILL_FOCUS, wxFocusEventHandler( GlobalconfDialog::updategc ), NULL, this );
	okbutton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( GlobalconfDialog::OK ), NULL, this );
	m_button6->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( GlobalconfDialog::cancel ), NULL, this );
}

GlobalconfDialog::~GlobalconfDialog()
{
	// Disconnect Events
	Gas->Disconnect( wxEVT_CHAR, wxKeyEventHandler( GlobalconfDialog::updategc ), NULL, this );
	Gas->Disconnect( wxEVT_KEY_DOWN, wxKeyEventHandler( GlobalconfDialog::updategc ), NULL, this );
	p0->Disconnect( wxEVT_KILL_FOCUS, wxFocusEventHandler( GlobalconfDialog::updategc ), NULL, this );
	p0->Disconnect( wxEVT_LEAVE_WINDOW, wxMouseEventHandler( GlobalconfDialog::updategc ), NULL, this );
	T0->Disconnect( wxEVT_KILL_FOCUS, wxFocusEventHandler( GlobalconfDialog::updategc ), NULL, this );
	freq->Disconnect( wxEVT_KILL_FOCUS, wxFocusEventHandler( GlobalconfDialog::updategc ), NULL, this );
	okbutton->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( GlobalconfDialog::OK ), NULL, this );
	m_button6->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( GlobalconfDialog::cancel ), NULL, this );
	
}

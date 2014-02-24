#include "gasvar.h"

namespace gases {
var gasvar::rho(const var& T,const var& p) {
	var result(T);
	vd rhotd=rho(T.tdata(),p.tdata());
	result.settdata(rhotd);
	return result;
}
var gasvar::p(const var& T,const var& rho) {
	var result(T);
	vd rhotdata=rho.tdata();
	vd Ttdata=T.tdata();
	vd restd=p(Ttdata,rhotdata);
	result.settdata(restd);
	return result;
}
var gasvar::cp(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=cp(tdata);
	result.settdata(restd);
	return result;
}
var gasvar::pr(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=pr(tdata);
	result.settdata(restd);
	return result;
}
var gasvar::h(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=h(tdata);
	result.settdata(restd);
	return result;
}
var gasvar::cv(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=cv(tdata);
	result.settdata(restd);
	return result;
}
var gasvar::e(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=e(tdata);
	result.settdata(restd);
	return result;
}
var gasvar::beta(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=beta(tdata);
	result.settdata(restd);
	return result;
}
var gasvar::gamma(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=gamma(tdata);
	result.settdata(restd);
	return result;
}
var gasvar::cm(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=cm(tdata);
	result.settdata(restd);
	return result;
}
var gasvar::mu(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=mu(tdata);
	result.settdata(restd);
	return result;
}
var gasvar::kappa(const var& T) {
	var result(T);
	vd tdata=T.tdata();
	vd restd=kappa(tdata);
	result.settdata(restd);
	return result;
}
}//Namespace gases

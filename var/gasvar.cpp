#include "gasvar.h"

namespace gases {

  var Gasvar::rho(const var& T,const var& p) {
    var result(T);
    vd Ttdata=T.tdata();
    vd ptdata=p.tdata();
    vd rhotd=Gas::rho(Ttdata,ptdata);
    result.settdata(rhotd);
    return result;
  }
  var Gasvar::p(const var& T,const var& rho) {
    var result(T);
    vd rhotdata=rho.tdata();
    vd Ttdata=T.tdata();
    vd restd=Gas::p(Ttdata,rhotdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::cp(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::cp(tdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::pr(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::pr(tdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::h(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::h(tdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::cv(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::cv(tdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::e(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::e(tdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::beta(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::beta(tdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::gamma(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::gamma(tdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::cm(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::cm(tdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::mu(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::mu(tdata);
    result.settdata(restd);
    return result;
  }
  var Gasvar::kappa(const var& T) {
    var result(T);
    vd tdata=T.tdata();
    vd restd=Gas::kappa(tdata);
    result.settdata(restd);
    return result;
  }
}//Namespace gases

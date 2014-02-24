//#include "lintube.h"
#include <boost/python.hpp>
#include <iostream>
using std::cout;
using std::endl;
char const* greet() {
	//cout << pi << endl;
	//cout << pi << endl;
    return "hello, world";
}


BOOST_PYTHON_MODULE(lintube_boost) {
    using namespace boost::python;
    def("greet", greet);
}

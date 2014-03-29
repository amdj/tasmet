#include <Python.h>
// First test of an implementation of a tube



static PyObject * Jacobian(PyObject *self, PyObject *args)
{
  return NULL;
}
static PyObject* Error(PyObject* self,PyObject *args){
  
  return NULL;
}


static PyMethodDef TubeMethods[] = {
  {"Jac",  Jacobian, METH_VARARGS,"Compute the Jacobian"},
  {"Err",  Error   , METH_VARARGS,"Compute the Error"}
  // {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initpytube(void)
{
  PyObject *m;

  m = Py_InitModule("pytube", TubeMethods);
  if (m == NULL)
    return;

}

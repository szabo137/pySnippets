/* Example of wrapping cos function from math.h with the python-C-API*/

#include <Python.h>
#include <math.h>

/* wrapped cos function */
static PyObject* cos_func(PyObject* self, PyObject* args)
{
	double value;
	double answer;

	//parse the input from py float to c double
	if (!PyArg_ParseTuple(args, "d", &value))
		return NULL;
	/*if the above function returns -1, an appropriate python exception will have been set and the function simply returns NULL*/

	//call cos from math.h
	answer = cos(cos(cos(cos(value))));

	//parse the output from c double to python float
	return Py_BuildValue("f",answer);
}

/* define function in module */

static PyMethodDef CosMethods[] =
{
	{"cos_func",cos_func,METH_VARARGS,"Evaluates the cos"},
	{NULL,NULL,0,NULL}
};

/* only for python 2 (for python 3, see docs) */
PyMODINIT_FUNC
initcos_module(void)
{
	(void) Py_InitModule("cos_module",CosMethods);
}

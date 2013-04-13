#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject* diffsim_calc_simulation(PyObject* self, PyObject* args);
PyMODINIT_FUNC PyInit_diffsim(void);

static PyMethodDef diffsimMethods[] = {
    {"calcsimulation", diffsim_calc_simulation, METH_VARARGS,
        "Executes the diffusion simulation wtih a set of parameters"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef diffsimmodule = {
    PyModuleDef_HEAD_INIT,
    "diffsim",
    "Accelerated module for diffusion simulation",
    -1,
    diffsimMethods
};

static PyObject* SimulationUnstableError;

PyMODINIT_FUNC PyInit_diffsim()
{
    PyObject* m = PyModule_Create(&diffsimmodule);
    if (m == NULL)
        return NULL;

    SimulationUnstableError = PyErr_NewException("diffsim.SimulationUnstableError", NULL, NULL);
    Py_INCREF(SimulationUnstableError);
    PyModule_AddObject(m, "SimulationUnstableError", SimulationUnstableError);

    import_array();

    return m;
}

static PyObject* diffsim_calc_simulation(PyObject* self, PyObject* args)
{
    PyObject* DVectorObj;
    PyObject* RVectorObj;
    int nIV;
    PyObject* initCondObj;
    int ndt, ndx;
    int dt, dx;
    double r;

    if (!PyArg_ParseTuple(args, "OOiOiiiid", &DVectorObj, &RVectorObj,
                                                &nIV, &initCondObj,
                                                &ndt, &ndx, &dt, &dx, &r))
    {
        return NULL;
    }
    if (DVectorObj == NULL || RVectorObj == NULL || initCondObj == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "One of the required array objects was null");
        return NULL;
    }

    PyObject* DVectorNewObj = PyArray_FROM_OTF(DVectorObj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject* RVectorNewObj = PyArray_FROM_OTF(RVectorObj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject* initCondNewObj = PyArray_FROM_OTF(initCondObj, NPY_DOUBLE, NPY_IN_ARRAY);

    if (DVectorObj == NULL || RVectorObj == NULL || initCondObj == NULL)
    {
        Py_XDECREF(DVectorNewObj);
        Py_XDECREF(RVectorNewObj);
        Py_XDECREF(initCondNewObj);

        PyErr_SetString(PyExc_ValueError, "One of the required array objects was null");
        return NULL;
    }

    double* DVectorC = PyArray_DATA(DVectorNewObj);
    double* RVectorC = PyArray_DATA(RVectorNewObj);
    double* initCondC = PyArray_DATA(initCondNewObj);
    //We just make the wonderful assumption that their len == nIV; don't fuck this precondition up!

    npy_intp dims[1];
    dims[0] = ndx;
    PyObject* returnObj = PyArray_Simple_New(1, dims, NPY_DOUBLE);


    Py_DECREF(DVectorNewObj);
    Py_DECREF(RVectorNewObj);
    Py_DECREF(initCondNewObj);
}

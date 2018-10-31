/*
 * Copyright (C) 2016  Leo Singer
 *
 * Fake implementations of some parts of the Python C API functions to aid
 * 2 to 3 porting. Named after the 'six' package that serves a similar
 * purpose for pure Python code (https://pypi.python.org/pypi/six).
 *
 *
 * Copied from https://github.com/lscsoft/lalsuite/blob/master/gnuscripts/six.h
 * but re-released under a 3-clause BSD style license - see LICENSE.rst
 */

#include <Python.h>

#ifndef _SIX_H
#define _SIX_H

#if PY_MAJOR_VERSION < 3
#ifdef PyMODINIT_FUNC
#undef PyMODINIT_FUNC
#define PyMODINIT_FUNC PyObject*
#endif

#define SIX_COMPAT_MODULE(name) _SIX_COMPAT_MODULE(name)
#define _SIX_COMPAT_MODULE(name) \
void init##name(void); /* Silence -Wmissing-prototypes */ \
void init##name(void) { \
    PyInit_##name(); \
}

typedef struct PyModuleDef {
    int m_base;
    const char *m_name;
    const char *m_doc;
    Py_ssize_t m_size;
    PyMethodDef *m_methods;
    inquiry m_reload;
    traverseproc m_traverse;
    inquiry m_clear;
    freefunc m_free;
} PyModuleDef;

#define PyModuleDef_HEAD_INIT 0

static PyObject *PyModule_Create(PyModuleDef *def)
{
    (void)PyModule_Create; /* Suppress unused function warning */

    if (def->m_size != -1)
    {
        PyErr_SetString(PyExc_NotImplementedError,
            "Python 2 does not support per-module state.");
        return NULL;
    }
    return Py_InitModule3(def->m_name, def->m_methods, def->m_doc);
}

static PyObject *PyModule_Create2(PyModuleDef *def, int module_api_version)
{
    (void)PyModule_Create2; /* Suppress unused function warning */

    if (def->m_size != -1)
    {
        PyErr_SetString(PyExc_NotImplementedError,
            "Python 2 does not support per-module state.");
        return NULL;
    }
    return Py_InitModule4(def->m_name, def->m_methods, def->m_doc,
        NULL, module_api_version);
}

#ifdef NUMPY_IMPORT_ARRAY_RETVAL
#undef NUMPY_IMPORT_ARRAY_RETVAL
#endif
#define NUMPY_IMPORT_ARRAY_RETVAL NULL

#ifdef NUMPY_IMPORT_UMATH_RETVAL
#undef NUMPY_IMPORT_UMATH_RETVAL
#endif
#define NUMPY_IMPORT_UMATH_RETVAL NULL

#else /* PY_MAJOR_VERSION >= 3 */
#define SIX_COMPAT_MODULE(name)
#endif

#endif /* _SIX_H */

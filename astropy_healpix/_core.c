#include <Python.h>
#include <math.h>
#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>
#include "healpix.h"
#include "healpix-utils.h"
#include "interpolation.h"

/* FIXME: We need npy_set_floatstatus_invalid(), but unlike most of the Numpy
 * C API it is only available on some platforms if you explicitly link against
 * Numpy, which is not typically done for building C extensions. This bundled
 * header file contains a static definition of _npy_set_floatstatus_invalid().
 */
#include "ieee754.h"

/* FIXME:
 * The Numpy C-API defines PyArrayDescr_Type as:
 *
 *   #define PyArrayDescr_Type (*(PyTypeObject *)PyArray_API[3])
 *
 * and then in some places we need to take its address, &PyArrayDescr_Type.
 * This is fine in GCC 10 and Clang, but earlier versions of GCC complain:
 *
 *   error: dereferencing pointer to incomplete type 'PyTypeObject'
 *   {aka 'struct _typeobject'}
 *
 * As a workaround, provide a faux forward declaration for PyTypeObject.
 * See https://github.com/numpy/numpy/issues/16970.
 *
 * Drop this when supporting gcc < 10 becomes irrelevant.
 */
#ifndef PYPY_VERSION
struct _typeobject {
    int _placeholder;
};
#endif


#define INVALID_INDEX (-1)

/* Data structure for storing function pointers for routines that are specific
 * to the HEALPix ordering scheme. When we create the ufuncs using
 * PyUFunc_FromFuncAndData, we will set them up to pass a pointer to this
 * data structure through as the void *data parameter for the loop functions
 * defined below. */
typedef struct {
    int64_t (*order_to_xy)(int64_t, int);
    int64_t (*xy_to_order)(int64_t, int);
} order_funcs;

static order_funcs
    nested_order_funcs = {healpixl_nested_to_xy, healpixl_xy_to_nested},
    ring_order_funcs   = {healpixl_ring_to_xy, healpixl_xy_to_ring};

static void *no_ufunc_data[]     = {NULL},
            *nested_ufunc_data[] = {&nested_order_funcs},
            *ring_ufunc_data[]   = {&ring_order_funcs};


static int64_t nside2npix(int nside)
{
    return 12 * ((int64_t) nside) * ((int64_t) nside);
}


static int pixel_nside_valid(int64_t pixel, int nside)
{
    return pixel >= 0 && pixel < nside2npix(nside);
}


static void healpix_to_lonlat_loop(
    char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
{
    order_funcs *funcs = data;
    npy_intp i, n = dimensions[0];

    for (i = 0; i < n; i ++)
    {
        int64_t pixel = *(int64_t *) &args[0][i * steps[0]];
        int     nside = *(int *)     &args[1][i * steps[1]];
        double  dx    = *(double *)  &args[2][i * steps[2]];
        double  dy    = *(double *)  &args[3][i * steps[3]];
        double  *lon  =  (double *)  &args[4][i * steps[4]];
        double  *lat  =  (double *)  &args[5][i * steps[5]];
        int64_t xy = INVALID_INDEX;

        if (pixel_nside_valid(pixel, nside))
            xy = funcs->order_to_xy(pixel, nside);

        if (xy >= 0)
            healpixl_to_radec(xy, nside, dx, dy, lon, lat);
        else
        {
            *lon = *lat = NPY_NAN;
            _npy_set_floatstatus_invalid();
        }
    }
}


static void lonlat_to_healpix_loop(
    char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
{
    order_funcs *funcs = data;
    npy_intp i, n = dimensions[0];

    for (i = 0; i < n; i ++)
    {
        double  lon    = *(double *)  &args[0][i * steps[0]];
        double  lat    = *(double *)  &args[1][i * steps[1]];
        int     nside  = *(int *)     &args[2][i * steps[2]];
        int64_t *pixel =  (int64_t *) &args[3][i * steps[3]];
        double  *dx    =  (double *)  &args[4][i * steps[4]];
        double  *dy    =  (double *)  &args[5][i * steps[5]];
        int64_t xy = INVALID_INDEX;

        if (isfinite(lon) && isfinite(lat))
            xy = radec_to_healpixlf(lon, lat, nside, dx, dy);
        if (xy >= 0)
            *pixel = funcs->xy_to_order(xy, nside);
        else {
            *pixel = INVALID_INDEX;
            *dx = *dy = NPY_NAN;
            _npy_set_floatstatus_invalid();
        }
    }
}


static void healpix_to_xyz_loop(
    char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
{
    order_funcs *funcs = data;
    npy_intp i, n = dimensions[0];

    for (i = 0; i < n; i ++)
    {
        int64_t pixel = *(int64_t *) &args[0][i * steps[0]];
        int     nside = *(int *)     &args[1][i * steps[1]];
        double  dx    = *(double *)  &args[2][i * steps[2]];
        double  dy    = *(double *)  &args[3][i * steps[3]];
        double  *x    =  (double *)  &args[4][i * steps[4]];
        double  *y    =  (double *)  &args[5][i * steps[5]];
        double  *z    =  (double *)  &args[6][i * steps[6]];
        int64_t xy = INVALID_INDEX;

        if (pixel_nside_valid(pixel, nside))
            xy = funcs->order_to_xy(pixel, nside);

        if (xy >= 0)
            healpixl_to_xyz(xy, nside, dx, dy, x, y, z);
        else
        {
            *x = *y = *z = NPY_NAN;
            _npy_set_floatstatus_invalid();
        }
    }
}


static void xyz_to_healpix_loop(
    char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
{
    order_funcs *funcs = data;
    npy_intp i, n = dimensions[0];

    for (i = 0; i < n; i ++)
    {
        double  x      = *(double *)  &args[0][i * steps[0]];
        double  y      = *(double *)  &args[1][i * steps[1]];
        double  z      = *(double *)  &args[2][i * steps[2]];
        int     nside  = *(int *)     &args[3][i * steps[3]];
        int64_t *pixel =  (int64_t *) &args[4][i * steps[4]];
        double  *dx    =  (double *)  &args[5][i * steps[5]];
        double  *dy    =  (double *)  &args[6][i * steps[6]];
        int64_t xy = INVALID_INDEX;

        if (isfinite(x) && isfinite(y) && isfinite(z)) {
            /* xyztohealpixlf expects a unit vector */
            double norm = sqrt(x*x + y*y + z*z);
            x /= norm;
            y /= norm;
            z /= norm;

            xy = xyztohealpixlf(x, y, z, nside, dx, dy);
        }
        if (xy >= 0)
            *pixel = funcs->xy_to_order(xy, nside);
        else {
            *pixel = INVALID_INDEX;
            *dx = *dy = NPY_NAN;
            _npy_set_floatstatus_invalid();
        }
    }
}


static void nested_to_ring_loop(
    char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
{
    npy_intp i, n = dimensions[0];

    for (i = 0; i < n; i ++)
    {
        int64_t nested = *(int64_t *) &args[0][i * steps[0]];
        int     nside  = *(int *)     &args[1][i * steps[1]];
        int64_t *ring  =  (int64_t *) &args[2][i * steps[2]];
        int64_t xy = INVALID_INDEX;

        if (pixel_nside_valid(nested, nside))
            xy = healpixl_nested_to_xy(nested, nside);
        if (xy >= 0)
            *ring = healpixl_xy_to_ring(xy, nside);
        else {
            *ring = INVALID_INDEX;
            _npy_set_floatstatus_invalid();
        }
    }
}


static void ring_to_nested_loop(
    char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
{
    npy_intp i, n = dimensions[0];

    for (i = 0; i < n; i ++)
    {
        int64_t ring   = *(int64_t *) &args[0][i * steps[0]];
        int     nside  = *(int *)     &args[1][i * steps[1]];
        int64_t *nested = (int64_t *) &args[2][i * steps[2]];
        int64_t xy = INVALID_INDEX;

        if (pixel_nside_valid(ring, nside))
            xy = healpixl_ring_to_xy(ring, nside);
        if (xy >= 0)
            *nested = healpixl_xy_to_nested(xy, nside);
        else {
            *nested = INVALID_INDEX;
            _npy_set_floatstatus_invalid();
        }
    }
}


static void bilinear_interpolation_weights_loop(
    char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
{
    npy_intp i, n = dimensions[0];

    for (i = 0; i < n; i ++)
    {
        double  lon    = *(double *)  &args[0][i * steps[0]];
        double  lat    = *(double *)  &args[1][i * steps[1]];
        int     nside  = *(int *)     &args[2][i * steps[2]];
        int64_t indices[4];
        double weights[4];
        int j;

        interpolate_weights(lon, lat, indices, weights, nside);
        for (j = 0; j < 4; j ++)
        {
            *(int64_t *) &args[3 + j][i * steps[3 + j]] = indices[j];
            *(double *)  &args[7 + j][i * steps[7 + j]] = weights[j];
        }
    }
}


static void neighbours_loop(
    char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
{
    order_funcs *funcs = data;
    npy_intp i, n = dimensions[0];

    for (i = 0; i < n; i ++)
    {
        int64_t pixel = *(int64_t *) &args[0][i * steps[0]];
        int     nside = *(int *)     &args[1][i * steps[1]];
        int64_t neighbours[] = {
            INVALID_INDEX, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX,
            INVALID_INDEX, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
        int j;
        int64_t xy = INVALID_INDEX;

        if (pixel_nside_valid(pixel, nside))
            xy = funcs->order_to_xy(pixel, nside);
        if (xy >= 0)
            healpixl_get_neighbours(xy, neighbours, nside);

        for (j = 0; j < 8; j ++)
        {
            int k = 4 - j;
            if (k < 0)
                k += 8;
            xy = neighbours[k];
            if (xy >= 0)
                pixel = funcs->xy_to_order(xy, nside);
            else {
                pixel = INVALID_INDEX;
                _npy_set_floatstatus_invalid();
            }
            *(int64_t *) &args[2 + j][i * steps[2 + j]] = pixel;
        }
    }
}


static PyObject *healpix_cone_search(
    PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *result;
    static char *kws[] = {"lon", "lat", "radius", "nside", "order", NULL};
    double lon, lat, radius;
    int nside;
    char *order;
    int64_t *indices, n_indices;
    int64_t *result_data;
    npy_intp dims[1];

    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "dddis", kws, &lon, &lat, &radius, &nside, &order))
        return NULL;

    n_indices = healpix_rangesearch_radec_simple(
        lon, lat, radius, nside, 0, &indices);
    if (!indices)
    {
        PyErr_SetString(
            PyExc_RuntimeError, "healpix_rangesearch_radec_simple failed");
        return NULL;
    }

    dims[0] = n_indices;
    result = PyArray_SimpleNew(1, dims, NPY_INT64);
    if (result)
    {
        result_data = PyArray_DATA((PyArrayObject *) result);

        if (strcmp(order, "nested") == 0)
        {
            int i;
            for (i = 0; i < n_indices; i ++)
                result_data[i] = healpixl_xy_to_nested(indices[i], nside);
        } else {
            int i;
            for (i = 0; i < n_indices; i ++)
                result_data[i] = healpixl_xy_to_ring(indices[i], nside);
        }
    }

    free(indices);
    return result;
}


static PyMethodDef methods[] = {
    {"healpix_cone_search", (PyCFunction) healpix_cone_search, METH_VARARGS | METH_KEYWORDS, NULL},
    {NULL, NULL, 0, NULL}
};


static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_core", NULL, -1, methods
};

static PyUFuncGenericFunction
    healpix_to_lonlat_loops             [] = {healpix_to_lonlat_loop},
    lonlat_to_healpix_loops             [] = {lonlat_to_healpix_loop},
    healpix_to_xyz_loops                [] = {healpix_to_xyz_loop},
    xyz_to_healpix_loops                [] = {xyz_to_healpix_loop},
    nested_to_ring_loops                [] = {nested_to_ring_loop},
    ring_to_nested_loops                [] = {ring_to_nested_loop},
    bilinear_interpolation_weights_loops[] = {bilinear_interpolation_weights_loop},
    neighbours_loops                    [] = {neighbours_loop};

static char
    healpix_to_lonlat_types[] = {
        NPY_INT64, NPY_INT, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE},
    lonlat_to_healpix_types[] = {
        NPY_DOUBLE, NPY_DOUBLE, NPY_INT, NPY_INT64, NPY_DOUBLE, NPY_DOUBLE},
    healpix_to_xyz_types[] = {
        NPY_INT64, NPY_INT, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE},
    xyz_to_healpix_types[] = {
        NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_INT, NPY_INT64, NPY_DOUBLE, NPY_DOUBLE},
    healpix_to_healpix_types[] = {
        NPY_INT64, NPY_INT, NPY_INT64},
    bilinear_interpolation_weights_types[] = {
        NPY_DOUBLE, NPY_DOUBLE, NPY_INT,
        NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64,
        NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE},
    neighbours_types[] = {
        NPY_INT64, NPY_INT,
        NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64,
        NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64};


PyMODINIT_FUNC PyInit__core(void)
{
    PyObject *module;

    import_array();
    import_umath();

    module = PyModule_Create(&moduledef);

    PyModule_AddObject(
        module, "healpix_nested_to_lonlat", PyUFunc_FromFuncAndData(
            healpix_to_lonlat_loops, nested_ufunc_data,
            healpix_to_lonlat_types, 1, 4, 2, PyUFunc_None,
            "healpix_nested_to_lonlat", NULL, 0));

    PyModule_AddObject(
        module, "healpix_ring_to_lonlat", PyUFunc_FromFuncAndData(
            healpix_to_lonlat_loops, ring_ufunc_data,
            healpix_to_lonlat_types, 1, 4, 2, PyUFunc_None,
            "healpix_ring_to_lonlat", NULL, 0));

    PyModule_AddObject(
        module, "lonlat_to_healpix_nested", PyUFunc_FromFuncAndData(
            lonlat_to_healpix_loops, nested_ufunc_data,
            lonlat_to_healpix_types, 1, 3, 3, PyUFunc_None,
            "lonlat_to_healpix_nested", NULL, 0));

    PyModule_AddObject(
        module, "lonlat_to_healpix_ring", PyUFunc_FromFuncAndData(
            lonlat_to_healpix_loops, ring_ufunc_data,
            lonlat_to_healpix_types, 1, 3, 3, PyUFunc_None,
            "lonlat_to_healpix_ring", NULL, 0));

    PyModule_AddObject(
        module, "healpix_nested_to_xyz", PyUFunc_FromFuncAndData(
            healpix_to_xyz_loops, nested_ufunc_data,
            healpix_to_xyz_types, 1, 4, 3, PyUFunc_None,
            "healpix_nested_to_xyz", NULL, 0));

    PyModule_AddObject(
        module, "healpix_ring_to_xyz", PyUFunc_FromFuncAndData(
            healpix_to_xyz_loops, ring_ufunc_data,
            healpix_to_xyz_types, 1, 4, 3, PyUFunc_None,
            "healpix_ring_to_xyz", NULL, 0));

    PyModule_AddObject(
        module, "xyz_to_healpix_nested", PyUFunc_FromFuncAndData(
            xyz_to_healpix_loops, nested_ufunc_data,
            xyz_to_healpix_types, 1, 4, 3, PyUFunc_None,
            "xyz_to_healpix_nested", NULL, 0));

    PyModule_AddObject(
        module, "xyz_to_healpix_ring", PyUFunc_FromFuncAndData(
            xyz_to_healpix_loops, ring_ufunc_data,
            xyz_to_healpix_types, 1, 4, 3, PyUFunc_None,
            "xyz_to_healpix_ring", NULL, 0));

    PyModule_AddObject(
        module, "nested_to_ring", PyUFunc_FromFuncAndData(
            nested_to_ring_loops, no_ufunc_data,
            healpix_to_healpix_types, 1, 2, 1, PyUFunc_None,
            "nested_to_ring", NULL, 0));

    PyModule_AddObject(
        module, "ring_to_nested", PyUFunc_FromFuncAndData(
            ring_to_nested_loops, no_ufunc_data,
            healpix_to_healpix_types, 1, 2, 1, PyUFunc_None,
            "ring_to_nested", NULL, 0));

    PyModule_AddObject(
        module, "bilinear_interpolation_weights", PyUFunc_FromFuncAndData(
            bilinear_interpolation_weights_loops, no_ufunc_data,
            bilinear_interpolation_weights_types, 1, 3, 8, PyUFunc_None,
            "bilinear_interpolation_weights", NULL, 0));

    PyModule_AddObject(
        module, "neighbours_nested", PyUFunc_FromFuncAndData(
            neighbours_loops, nested_ufunc_data,
            neighbours_types, 1, 2, 8, PyUFunc_None,
            "neighbours_nested", NULL, 0));

    PyModule_AddObject(
        module, "neighbours_ring", PyUFunc_FromFuncAndData(
            neighbours_loops, ring_ufunc_data,
            neighbours_types, 1, 2, 8, PyUFunc_None,
            "neighbours_ring", NULL, 0));

    return module;
}

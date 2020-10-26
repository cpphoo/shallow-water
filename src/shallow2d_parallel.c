#include <string.h>
#include <math.h>

//ldoc on
/**
 * ## Implementation
 *
 * The actually work of computing the fluxes and speeds is done
 * by local (`static`) helper functions that take as arguments
 * pointers to all the individual fields.  This is helpful to the
 * compilers, since by specifying the `restrict` keyword, we are
 * promising that we will not access the field data through the
 * wrong pointer.  This lets the compiler do a better job with
 * vectorization.
 */


static const float g = 9.8;                      // The gravitational constant.


static
void shallow2dv_flux(float* restrict FU_h,
                     float* restrict FU_hu,
                     float* restrict FU_hv,
                     float* restrict GU_h,
                     float* restrict GU_hu,
                     float* restrict GU_hv,
                     const float* restrict h,
                     const float* restrict hu,
                     const float* restrict hv,
                     float g,
                     int ncell)
{
  /* Description: This function calculates the components of FU and GU (the
  first 6 arguments) using the components of U (the next three arguments) and g.
  In particular, we use the fist ncell elements of h, hu, and hv to calculate
  the first ncell elements of FU_h, FU_hu, FU_hv, and GU_h, GU_hu, GU_hv.

  Arguments:
  FU_h, FU_hu, FU_hv - the three componnets of FU.

  GU_h, GU_hu, GU_hv - the three components of GU.

  h, hu, hv - the three components of U.

  g - the gravitational constant

  ncell - the number of cells that you want to operate on. */

  // If we look back at the equations on page 3, we can see that the h
  // component of FU is hu, and the h component of GU is hv. Thus, we can
  // just copy the hu and hv into FU_h and GU_h, respectively.
  memcpy(FU_h, hu, ncell * sizeof(float));
  memcpy(GU_h, hv, ncell * sizeof(float));

  // Calculate FU_hu, FU_hv, GU_hu, and GU_hv.
  for (int i = 0; i < ncell; ++i) {
    // first find h, hu and hv for the current cell. Use h to find 1/h (inv_h)
    float hi = h[i], hui = hu[i], hvi = hv[i];
    float inv_h = 1.0f/hi;

    // Now uses these to calculate FU_hu, FU_hv, GU_hu, and GU_hv for this cell.
    // This simply implements the equations on page 3 of the document.
    FU_hu[i] = hui*hui*inv_h + (0.5f*g)*hi*hi;
    FU_hv[i] = hui*hvi*inv_h;
    GU_hu[i] = hui*hvi*inv_h;
    GU_hv[i] = hvi*hvi*inv_h + (0.5f*g)*hi*hi;
  } // for (int i = 0; i < ncell; ++i) {
} // void shallow2dv_flux(float* restrict FU_h,...


static
void shallow2dv_speed(float* restrict cxy,
                      const float* restrict h,
                      const float* restrict hu,
                      const float* restrict hv,
                      float g,
                      int ncell)
{
  /* Description: This function calculates the maximum velocity in the x (cx)
  and y (cy) directions. We start off with an initial value for cx and cy
  (which are passed to the function through the cxy argumnet). This function
  calculates
    cx = max{ cx, cx_1, cx_2,... cx_ncell}
    cy = max{ cy, cy_1,... cy_ncell}
  Where for each i in {1,2... ncell}, cxi denotes the wave velocity in the x
  direction in the ith cell and cyi denotes the wave velocity in the y direction
  in the ith cell.

  What are the arguments?
  cxy - a 2 element array. On input, this holds some initial values. After the
  function runs, cxy[0] will hold max{cxy[0], cx_1,... cx_ncell} (see above)

  h, hu, hv - the components of U.

  g - the gravitational constant

  ncell - the number of cells in h, hu, hv that we want to operate on. */

  // set up local variables to keep track of cx and cy .
  float cx = cxy[0];
  float cy = cxy[1];

  #pragma omp parallel for
  // Cycle through the first ncell elements of h, hu, and hv.
  for (int i = 0; i < ncell; ++i) {
    // find h, 1/h (h_inv) and sqrt(gh) for this cell.
    float hi = h[i];
    float inv_hi = 1.0f/h[i];
    float root_gh = sqrtf(g * hi);

    // cxi = |u[i]| + sqrt(gh).
    // cyi = |v[i]| + sqrt(gh).
    float cxi = fabsf(hu[i] * inv_hi) + root_gh;
    float cyi = fabsf(hv[i] * inv_hi) + root_gh;

    // Check if velocities in this cell exceded the running maxumum (cx and cy)
    if (cx < cxi) cx = cxi;
    if (cy < cyi) cy = cyi;
  } // for (int i = 0; i < ncell; ++i) {

  // store the final (maximum) cx and cy back in the cxy argument.
  cxy[0] = cx;
  cxy[1] = cy;
} // void shallow2dv_speed(float* restrict cxy,...



void shallow2d_flux(float* FU,
                    float* GU,
                    const float* U,
                    int ncell,
                    int field_stride)
{
  /* Description: This function is a a wrapper for shallow2dv_flux. All it does
  is call shallow2dv_flux.

  What are the arguments?
  FU, GU, U - field_stride by 3 arrays stored in column major order. Each row
  holds the three components of U for a particular cell (see page 3 of the
  document). The first holds the "h" component, the second holds the "hu"
  component and the third holds the "hv" component.

  ncell - the number of cells we want to operate on.

  field_stride - the distance in memory between successive fields (components) of
  FU, GU, and U. For example, the distance in memory between the h and hu
  componenets of U for the ith cell is field_stride. */

  shallow2dv_flux(FU,
                  FU + field_stride,
                  FU + 2*field_stride,
                  GU,
                  GU + field_stride,
                  GU + 2*field_stride,
                  U,
                  U + field_stride,
                  U + 2*field_stride,
                  g,
                  ncell);
} // void shallow2d_flux(float* FU,...



void shallow2d_speed( float* cxy,
                      const float* U,
                      int ncell,
                      int field_stride)
{
  /* Description: This function is a wrapper for shallow2dv_speed. All it does
  is call shallow2dv_speed.

  What are the arguments?
  cxy - a 2 element array. On input, cxy[0] holds some initial value. On output,
  cxy[0] holds the maximum of the initial value of cxy[0] and the value of the
  x velocity in the first ncell cells of U. cxy[1] does the same for the y
  velocity.

  U - a field cell by c array stored in column major order. Each row holds the
  three components of U for a particular cell (see page 3 of the document).

  ncell - the number of cells we want to operate on .

  field_stride - the distance in memory between successive fields (components)
  of U. For example, the distance in memory between the h and hu
  componenets of U for the ith cell is field_stride. */

  shallow2dv_speed( cxy,
                    U,
                    U + field_stride,
                    U + 2*field_stride,
                    g,
                    ncell);
} // void shallow2d_speed(float* cxy,...

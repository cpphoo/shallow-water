#include "stepper_parallel.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <omp.h>
#include <stdio.h>

//ldoc on
/**
 * ## Implementation
 *
 * ### Structure allocation
 */

central2d_t* central2d_init(float grid_width,
                            float grid_height,
                            int nx,
                            int ny,
                            int nfield,
                            flux_t flux,
                            speed_t speed,
                            float cfl)
{
  /* Description: This function is essentially a constructor for a central2d
  structure. The parameters that it accepts are used to set up a central2d
  structure. This function returns a pointer to that structure.

  What are the arguments?
  grid_width - the length of the grid (of canonical cells) in the x direction.

  grid_height - the length of the grid (of canonical cells) in the y direction.

  nx - number of cells in the x direction.

  ny - number of cells in the y direction

  nfield - The of quantities that we keep track of in each cell. Equivalently,
  this is the number of components in FU, GU, or U for a given cell (3 in our
  case). Equivalently, this is the number of sub-arrays in FU, Gu, and U.

  flux - a pointer to a function to update F and G using U (see shallow2d.c)

  speed - pointer to a function that will calculate the maximum wave speed in
  the x and y directions (see shallow2d.c)

  cfl - dictates the the allowed time step. */

  // We extend to a four cell buffer to avoid BC comm on odd time steps
  int ng = 4;

  // allocate central2d object.
  central2d_t* sim = (central2d_t*) malloc(sizeof(central2d_t));

  // populate some of its members.
  sim->nx = nx;
  sim->ny = ny;
  sim->ng = ng;
  sim->nfield = nfield;
  sim->dx = grid_width/nx;
  sim->dy = grid_height/ny;
  sim->flux = flux;
  sim->speed = speed;
  sim->cfl = cfl;

  // Calculate nx_all, ny_all... use these to calculate N (the size of
  // U, V, FU, or GU)
  int nx_all = nx + 2*ng;
  int ny_all = ny + 2*ng;
  int ncells = nx_all*ny_all;
  int N  = nfield*ncells;

  // Here we allocate a single large block of memory. We split this block
  // into five pieces. The first four go to U, U_half, FU and GU, respectively
  // Each of these pieces have size N.
  // The last piece, which has size 6*nx_all, goes to scratch (which is
  // essentially used to hold "scratch" calculations that are used in other
  // calculations but which we don't need to keep track of).
  sim->U  = (float*) malloc((4*N + 6*nx_all)* sizeof(float));
  sim->U_half  = sim->U +   N;
  sim->FU      = sim->U + 2*N;
  sim->GU      = sim->U + 3*N;
  sim->scratch = sim->U + 4*N;

  return sim;
} // central2d_t* central2d_init(float w,



void central2d_free(central2d_t* sim)
{
  /* Description: This function is a destructor for central2d structures.

  What is the argument?
  sim - a pointer to the central2d object that you want to destroy. */

  free(sim->U);
  free(sim);
} // void central2d_free(central2d_t* sim)



int central2d_offset( central2d_t* sim,
                      int k,
                      int ix,
                      int iy)
{
  /* Description: this function is used to access different elements of U.
  U is an n_fields by nx_all by ny_all array. We think of U as a sequence of
  three nx_all by ny_all subarrays, each of whcih is stored in ROW MAJOR order.

  The first subarray (the first nx_all*ny_all elements of U) store the value of
  "h" for each cell in the grid (including ghost cells)

  The second subarray (the next nx_all*ny_all elements of U) store the value of
  "hu" for each cell in the grid (including ghost cells).

  The third subarray (the final nx_all*ny_all elements of U) store the value of
  "hv" for each cell in the grid (including ghost cells),

  What are the arguments?
  sim - a pointer to the central2d structure whose U array you want to access.

  k - the field (component) number you want to access. For us,
    k = 0 means "h" componet
    k = 1 means "hu" component
    k = 2 means "hv" component

  ix, iy - the x, y indicies of a canonical cell (within the canonical grid) that
  we want to access. ix should be in the set {0, 1... nx-1} and iy should be
  in the set {0, 1... ny-1} */

  // How does this work?
  //
  // Suppose that we want to find the hu component of a particular cell in the
  // grid. We know that U + nx_all*ny_all is the first index of the subarray
  // of U whcih stores the "hu" component of every cell in the grid (including
  // ghost cells). The first and last ng rows and columns of the sub array
  // correspond to ghost cells. Thus, the column of the sub array corresponding
  // to the iyth “canonical” column is iy + ng. Likewise, the row of the sub
  // array corresponding to the ixth “canonical” cell is ng + ix. Therefore, the
  // location (in memory, within the kth sub array of C) of the value of the "hu"
  // component of the (ix, iy) canonical cell is (ng + iy)*nx_all + (ng + ix).
  // (remember, each sub array is in ROW MAJOR order)
  //
  // Combining this with the fact that the first entry of the kth sub array of
  // U is at k*nx_all*ny_all, we can conclude that the desired quantity is at
  // the following location:
  //
  //    U + k*nx_all*ny_all + (ng+iy)*nx_all + (ng+ix).
  //
  // Which is exactly what this function returns.

  // fetch nx, ny, ng, use them to calculate nx_all and ny_all
  int nx = sim->nx, ny = sim->ny, ng = sim->ng;
  int nx_all = nx + 2*ng;
  int ny_all = ny + 2*ng;

  // Return the desired quantitity (See "How does this work?" above)
  return (k*ny_all + (ng + iy))*nx_all + (ng + ix);
} // int central2d_offset( central2d_t* sim,..


/**
 * ### Boundary conditions
 *
 * In finite volume methods, boundary conditions are typically applied by
 * setting appropriate values in ghost cells.  For our framework, we will
 * apply periodic boundary conditions; that is, waves that exit one side
 * of the domain will enter from the other side.
 *
 * We apply the conditions by assuming that the cells with coordinates
 * `nghost <= ix <= nx+nghost` and `nghost <= iy <= ny+nghost` are
 * "canonical", and setting the values for all other cells `(ix,iy)`
 * to the corresponding canonical values `(ix+p*nx,iy+q*ny)` for some
 * integers `p` and `q`.
 */

static inline
void copy_subgrid(float* restrict dst,
                 const float* restrict src,
                 int nx,
                 int ny,
                 int stride_dst,
                 int stride_src)
{
  /* Description: This copies an nx by ny part of src to dest.
  Both src and dst are assumed to be in ROW MAJOR order.

  What are the arguments?
  dst - the destination of the copy

  src - the source of the copy

  nx - the number of cells in the x direction to be copied.

  ny - the number of cells in the y direction to be copied.

  stride_dst - the distance (in memory) between successive elements of a column
  of dst. Equivalently, the size (in memory) of a row of dst.

  stride_src - the distance (in memory) between successive elements of a column
  of src. Equivalently, the size (in memory) of a row of src. */

  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      dst[ix + iy*stride_dst] = src[ix + iy*stride_src];
    } // for (int ix = 0; ix < nx; ++ix) {
  } // for (int iy = 0; iy < ny; ++iy) {
} // void copy_subgrid(float* restrict dst,...


void central2d_periodic(float* restrict U_global,
                        int nx,
                        int ny,
                        int ng,
                        int nfield)
{
  /* Description: This function applies periodic boundary conditions to U.
  U is an nfield by nx_all by ny_all array. We think of U as a sequence of
  nfield sub-arrays, each of size nx_all by ny_all. Each sub array is stored
  in ROW MAJOR order.

  What are the arguments?
  U - the U array of a central2d structure

  nx - number of canonical cells in the x direction of the grid

  ny - the number of canonical cells in the y direction of the grid

  ng - number of layers of ghost cells (the first and las ng rows of U are
  ghost cells. Likewise, the first and last ng columns of U are ghost cells).

  nfield - the number of fields/componenets/subarrays in U.
  */

  // Stride and number per field.
  int nx_all = nx + 2*ng;
  int ny_all = ny + 2*ng;
  int s = nx_all;
  int field_stride = nx_all*ny_all;

  /* Offsets of left, right, top, and bottom data blocks and ghost blocks

  Increasing the row index increases the y value. Increasing the column index
  increases the x value. Thus, the top of the grid corresponds to the last row.
  Likewise, the right of the grid corresponds to the last column.

  l denotes the address of the bottom right (min x and y indicies) of the cells
  in U which will become the left boundary of the ghost cells in U.

  lg denotes the address of the bottom right (min x and y indicies) of the left
  boundary cells in U.

  r, b, t, br, bg and tg are similar. If we think about it, the locations of
  l, r, b, t, lg, lr, bg, and tg should make sense. (draw a picture, it helps!) */
  int l = nx,   lg = 0;
  int r = ng,   rg = nx + ng;
  int t = ng*s, tg = (ny + ng)*s;
  int b = ny*s, bg = 0;

  // Copy data into ghost cells on each side
  for (int k = 0; k < nfield; ++k) {
    // Get the address of the kth subarray of U.
    float* Uk_global = U_global + k*field_stride;

    copy_subgrid(Uk_global + lg, Uk_global + l, ng,     ny_all, s, s);
    copy_subgrid(Uk_global + rg, Uk_global + r, ng,     ny_all, s, s);
    copy_subgrid(Uk_global + tg, Uk_global + t, nx_all, ng,     s, s);
    copy_subgrid(Uk_global + bg, Uk_global + b, nx_all, ng,     s, s);
  } // for (int k = 0; k < nfield; ++k) {
} // void central2d_periodic(float* restrict U,


void central2d_local_BC(float* restrict U,
                        int nx,
                        int ny,
                        int ng,
                        int nfield,
                        float* restrict U_global,
                        int nx_global,
                        int ny_global,
                        int xlow_local,
                        int ylow_local)
{
  /* Description: This function applies boundary conditions to U a piece of
  the partition of the global U array.
  U is an nfield by nx_all by ny_all array. We think of U as a sequence of
  nfield sub-arrays, each of size nx_all by ny_all. Each sub array is stored
  in ROW MAJOR order.

  What are the arguments?
  U - the U array of a central2d structure for a piece of the global grid.

  nx - number of canonical cells in the x direction of the grid

  ny - the number of canonical cells in the y direction of the grid

  ng - number of layers of ghost cells (the first and las ng rows of U are
  ghost cells. Likewise, the first and last ng columns of U are ghost cells).

  nfield - the number of fields/componenets/subarrays in U.

  U_global - the U array for the global grid.

  nx_global - the number of canonical cells in the x direction of the global
  grid.

  ny_global - the number of canonical cells in the y direction of the global
  grid.

  xlow_local, ylow_local - these tell us how U fits into U local. In
  particular, the first row and column of canonical cells in U corresponds to
  the xlow_local row and ylow_local column of the global U. */

  // Stride and number per field for local U.
  int nx_all = nx + 2*ng;
  int ny_all = ny + 2*ng;
  int stride_local = nx_all;
  int field_stride = nx_all*ny_all;

  // Stride and number per field for global U.
  int nx_global_all = nx_global + 2*ng;
  int ny_global_all = ny_global + 2*ng;
  int stride_global = nx_global_all;
  int field_stride_global = nx_global_all*ny_global_all;

  /* Offsets of left, right, top, and bottom data blocks and ghost blocks

  Increasing the row index increases the y value. Increasing the column index
  increases the x value. Thus, the top of the grid corresponds to the last row.
  Likewise, the right of the grid corresponds to the last column.

  l denotes the address of the bottom right (min x and y indicies) of the cells
  in U_global which will become the left boundary of the ghost cells in U.

  lg denotes the address of the bottom right (min x and y indicies) of the left
  boundary cells in U.

  r, b, t, rg, bg, and tg are similar.

  To get these, we first need to find the address of the cell of U_global (for a
  particular field/component/subarray) which corresponds to the bottom left
  ghost cell of U.

  To do this, let's think about which row and column this corresponds to. There
  are nx_global_all columns in the global grid (with ghost cells). The first ng
  of those columns hold ghost cells (which we used to apply BCs!). By
  definition, xlow_local tells us which column (of canoncal cells) in U_global
  corresponds to the first column of canonical cells of U. Thus, the first
  column of caonical cells in U coresponds to column xlow_local + ng in U_global
  The first column of ghost cells is ng columns to the left of that. Thus, the
  column index in U_global corresonding to the leftmost column of ghost cells in
  U is xlow_local + ng - ng = xlow_local.

  Using analagous logic, we can conclude that the bottom row of ghost cells in U
  corresponds to row ylow_local + ng - ng = ylow_local in U_global.

  Since there are nx_global_all entries per row of U_global (which is stored in
  ROW MAJOR order), the address in U_global corresponding to the bottom
  left most ghost cell in U is xlow_local + ylow_local*nx_global_all. We call
  this quantity offset. */
  int offset = xlow_local + ylow_local*nx_global_all;

  int l = offset;
  int lg = 0;

  int r = offset + (nx + ng);
  int rg = nx + ng;

  int t = offset + (ny + ng)*nx_global_all;
  int tg = (ny + ng)*nx_all;

  int b = offset;
  int bg = 0;

  for (int k = 0; k < nfield; ++k) {
    // Get the address of the kth subarray of U and U_global.
    float* Uk_local = U + k*field_stride;
    float* Uk_global = U_global + k*field_stride_global;

    copy_subgrid(Uk_local + lg, Uk_global + l, ng,     ny_all, stride_local, stride_global);
    copy_subgrid(Uk_local + rg, Uk_global + r, ng,     ny_all, stride_local, stride_global);
    copy_subgrid(Uk_local + tg, Uk_global + t, nx_all, ng,     stride_local, stride_global);
    copy_subgrid(Uk_local + bg, Uk_global + b, nx_all, ng,     stride_local, stride_global);
  } // for (int k = 0; k < nfield; ++k) {
} // void central2d_local_BC(float* restrict U,...



void central2d_local_to_global(float* restrict U,
                               int nx,
                               int ny,
                               int ng,
                               int nfield,
                               float* restrict U_global,
                               int nx_global,
                               int ny_global,
                               int xlow_local,
                               int ylow_local)
{
  /* Description: This function is used to map the first and last ng rows of
  canonical cells in U to their corresponding entries in U_global. These
  elements are the only elements that the threads communicate to one another.
  Thus, this function gets global U ready for the next time step.

  What are the arguments?
  U - the U array of a central2d structure for a piece of the global grid.

  nx - number of canonical cells in the x direction of the grid

  ny - the number of canonical cells in the y direction of the grid

  ng - number of layers of ghost cells (the first and las ng rows of U are
  ghost cells. Likewise, the first and last ng columns of U are ghost cells).

  nfield - the number of fields/componenets/subarrays in U.

  U_global - the U array for the global grid.

  nx_global - the number of canonical cells in the x direction of the global
  grid.

  ny_global - the number of canonical cells in the y direction of the global
  grid.

  xlow_local, ylow_local - these tell us how U fits into U local. In
  particular, the first row and column of canonical cells in U corresponds to
  the xlow_local row and ylow_local column of the global U. */

  // Stride and number per field for local U.
  int nx_all = nx + 2*ng;
  int ny_all = ny + 2*ng;
  int stride_local = nx_all;
  int field_stride = nx_all*ny_all;

  // Stried and number per field for global U.
  int nx_global_all = nx_global + 2*ng;
  int ny_global_all = ny_global + 2*ng;
  int stride_global = nx_global_all;
  int field_stride_global = nx_global_all*ny_global_all;

  /* Find the addresses of the top, bottom, left, and right most canonical
  cells within U and U_global.

  Increasing the row index increases the y value. Increasing the column index
  increases the x value. Thus, the top of the grid corresponds to the last row.
  Likewise, the right of the grid corresponds to the last column.

  l_local denotes the address of the bottom right (min x and y indicies) cell
  in U_global which corresponds to the bottom right canonical cell in U.

  l_global denotes the address of the bottom right caonical cell in U.

  r_local, r_global, b_local, b_global, t_local, and t_global are similar.

  To get these, we first need to find the address of the cell of U_global (for a
  particular field/component/subarray) which corresponds to the bottom left
  caonical cell in U local.

  To do this, let's think about which row and column this corresponds to. There
  are nx_global_all columns in the global grid (with ghost cells). The first ng
  of those columns hold ghost cells (which we used to apply BCs!). By
  definition, xlow_local tells us which column (of canoncal cells) in U_global
  corresponds to the first column of canonical cells of U. Thus, the first
  column of caonical cells in U coresponds to column xlow_local + ng in U_global

  Using analagous logic, we can conclude that the bottom row of caonical in U
  corresponds to row ylow_local + ng in U_global.

  Since there are nx_global_all entries per row of U_global (which is stored in
  ROW MAJOR order), the address in U_global corresponding to the bottom
  left most caonical cell in U is (xlow_local + ng) + (ylow_local + ng)*nx_global_all.
  We call this quantity offset_global.

  Similarly, we find the address within U of the bottom left most caonical cell.
  This quantity, which we call offset_local, is equal to ng + ng*nx_all (think
  about it) */
  int offset_global = (xlow_local + ng) + (ylow_local + ng)*nx_global_all;
  int offset_local = ng + ng*nx_all;

  int l_local  = offset_local;
  int l_global = offset_global;

  int r_local  = offset_local  + (nx - ng);
  int r_global = offset_global + (nx - ng);

  int t_local  = offset_local  + (ny - ng)*nx_all;
  int t_global = offset_global + (ny - ng)*nx_global_all;

  int b_local  = offset_local;
  int b_global = offset_global;

  for (int k = 0; k < nfield; ++k) {
    // Get the address of the kth subarray of U and U_global.
    float* Uk_local  = U + k*field_stride;
    float* Uk_global = U_global + k*field_stride_global;

    copy_subgrid(Uk_global + l_global, Uk_local + l_local, ng,     ny,     stride_global, stride_local);
    copy_subgrid(Uk_global + r_global, Uk_local + r_local, ng,     ny,     stride_global, stride_local);
    copy_subgrid(Uk_global + t_global, Uk_local + t_local, nx,     ng,     stride_global, stride_local);
    copy_subgrid(Uk_global + b_global, Uk_local + b_local, nx,     ng,     stride_global, stride_local);
  } // for (int k = 0; k < nfield; ++k) {
} // void central2d_local_to_global(float* restrict U,



void central2d_U_to_global_U(float* restrict U,
                             int nx,
                             int ny,
                             int ng,
                             int nfield,
                             float* restrict U_global,
                             int nx_global,
                             int ny_global,
                             int xlow_local,
                             int ylow_local)
{
  /* Description: This function moves U to global U.
  U is an nfield by nx_all by ny_all array. We think of U as a sequence of
  nfield sub-arrays, each of size nx_all by ny_all. Each sub array is stored
  in ROW MAJOR order.

  What are the arguments?
  U - the U array of a central2d structure for a piece of the global grid.

  nx - number of canonical cells in the x direction of the grid

  ny - the number of canonical cells in the y direction of the grid

  ng - number of layers of ghost cells (the first and las ng rows of U are
  ghost cells. Likewise, the first and last ng columns of U are ghost cells).

  nfield - the number of fields/componenets/subarrays in U.

  U_global - the U array for the global grid.

  nx_global - the number of canonical cells in the x direction of the global
  grid.

  ny_global - the number of canonical cells in the y direction of the global
  grid.

  xlow_local, ylow_local - these tell us how U fits into U local. In
  particular, the first row and column of canonical cells in U corresponds to
  the xlow_local row and ylow_local column of the global U. */

  // Stride and number per field for local U.
  int nx_all = nx + 2*ng;
  int ny_all = ny + 2*ng;
  int stride_local = nx_all;
  int field_stride = nx_all*ny_all;

  // Stried and number per field for global U.
  int nx_global_all = nx_global + 2*ng;
  int ny_global_all = ny_global + 2*ng;
  int stride_global = nx_global_all;
  int field_stride_global = nx_global_all*ny_global_all;

  /* Offsets of U in U_global.

  We need to find the address the cell of U_global (for a particular
  field/component/subarray) which corresponds to the bottom left caonical
  cell of U.

  To do this, let's think about which row and column this corresponds to. There
  are nx_global_all columns in the global grid (with ghost cells). The first ng
  of those columns hold ghost cells (which we used to apply BCs!). By
  definition, xlow_local tells us which column (of canoncal cells) in U_global
  corresponds to the first column of canonical cells of U. Thus, the first
  column of caonical cells in U coresponds to column xlow_local + ng in U_global

  Using analagous logic, we can conclude that the bottom row of canonical cells
  in U corresponds to row ylow_local + ng in U_global.

  Since there are nx_global_all entries per row of U_global (which is stored in
  ROW MAJOR order), the address in U_global corresponding to the bottom
  left most canonical cell in U is (xlow_local + ng) + (ylow_local + ng)*nx_global_all.
  We call this quantity offset_global.

  Similarly, we find the address within U of the bottom left most caonical cell.
  This quantity, which we call offset_local, is equal to ng + ng*nx_all (think
  about it) */
  int offset_global = (xlow_local + ng) + (ylow_local + ng)*nx_global_all;
  int offset_local = ng + ng*nx_all;

  for (int k = 0; k < nfield; ++k) {
    // Get the address of the kth subarray of U and U_global.
    float* Uk_local = U + k*field_stride;
    float* Uk_global = U_global + k*field_stride_global;

    // Copy the canonical cells of U to their corresponding locations in U_global
    copy_subgrid(Uk_global + offset_global, Uk_local + offset_local, nx, ny, stride_global, stride_local);
  } // for (int k = 0; k < nfield; ++k) {
} // void central2d_local_BC(float* restrict U,...



/**
 * ### Derivatives with limiters
 *
 * In order to advance the time step, we also need to estimate
 * derivatives of the fluxes and the solution values at each cell.
 * In order to maintain stability, we apply a limiter here.
 *
 * The minmod limiter *looks* like it should be expensive to computer,
 * since superficially it seems to require a number of branches.
 * We do something a little tricky, getting rid of the condition
 * on the sign of the arguments using the `copysign` instruction.
 * If the compiler does the "right" thing with `max` and `min`
 * for floating point arguments (translating them to branch-free
 * intrinsic operations), this implementation should be relatively fast.
 */


// Branch-free computation of minmod of two numbers times 2s
static inline
float xmin2s(float s, float a, float b) {
    float sa = copysignf(s, a);
    float sb = copysignf(s, b);
    float abs_a = fabsf(a);
    float abs_b = fabsf(b);
    float min_abs = (abs_a < abs_b ? abs_a : abs_b);
    return (sa+sb) * min_abs;
}


// Limited combined slope estimate
static inline
float limdiff(float um, float u0, float up) {
    const float theta = 2.0;
    const float quarter = 0.25;
    float du1 = u0 - um;   // Difference to left
    float du2 = up - u0;   // Difference to right
    float duc = up - um;   // Twice centered difference
    return xmin2s( quarter, xmin2s(theta, du1, du2), duc );
}


// Compute limited derivs
static inline
void limited_deriv1(float* restrict du,
                    const float* restrict u,
                    int ncell)
{
    for (int i = 0; i < ncell; ++i)
        du[i] = limdiff(u[i - 1], u[i], u[i + 1]);
}


// Compute limited derivs across stride
static inline
void limited_derivk(float* restrict du,
                    const float* restrict u,
                    int ncell, int stride)
{
    assert(stride > 0);
    for (int i = 0; i < ncell; ++i)
        du[i] = limdiff(u[i - stride], u[i], u[i + stride]);
}


/**
 * ### Advancing a time step
 *
 * Take one step of the numerical scheme.  This consists of two pieces:
 * a first-order corrector computed at a half time step, which is used
 * to obtain new $F$ and $G$ values; and a corrector step that computes
 * the solution at the full step.  For full details, we refer to the
 * [Jiang and Tadmor paper][jt].
 *
 * The `compute_step` function takes two arguments: the `io` flag
 * which is the time step modulo 2 (0 if even, 1 if odd); and the `dt`
 * flag, which actually determines the time step length.  We need
 * to know the even-vs-odd distinction because the Jiang-Tadmor
 * scheme alternates between a primary grid (on even steps) and a
 * staggered grid (on odd steps).  This means that the data at $(i,j)$
 * in an even step and the data at $(i,j)$ in an odd step represent
 * values at different locations in space, offset by half a space step
 * in each direction.  Every other step, we shift things back by one
 * mesh cell in each direction, essentially resetting to the primary
 * indexing scheme.
 *
 * We're slightly tricky in the corrector in that we write
 * $$
 *   v(i,j) = (s(i+1,j) + s(i,j)) - (d(i+1,j)-d(i,j))
 * $$
 * where $s(i,j)$ comprises the $u$ and $x$-derivative terms in the
 * update formula, and $d(i,j)$ the $y$-derivative terms.  This cuts
 * the arithmetic cost a little (not that it's that big to start).
 * It also makes it more obvious that we only need four rows worth
 * of scratch space.
 */


// Predictor half-step
static
void central2d_predict(float* restrict U_half,
                       float* restrict scratch,
                       const float* restrict U,
                       const float* restrict FU,
                       const float* restrict GU,
                       float dtcdx2,
                       float dtcdy2,
                       int nx_all,
                       int ny_all,
                       int nfield)
{
  /* Description: This function uses U, FU, and FG to calculate U_half (U on a
  grid which is staggered by a half step with respect to the grid that U is on).

  What are the arguments?
  U_half, scratch, U, FU, GU - the addresses of the corresponding members of
  a central2d structure. Each of U_half, U, FU, and GU are nfiled by nx_all
  by ny_all arrays (3d arrays!)

  dtcdx2 - (1/2)(dt/dx) (see central2d_step)

  dtcdy2 -(1/2)(dt/dy) (see central2d_step)

  nx_all - number of columns (x direction) of cells (including ghost cells).

  ny_all - number of rows (y direction) of cells (including ghost cells).

  nfield - the number of fields/subarrays/components in U, U_half, FU, FG. */


  float* restrict FUx = scratch;
  float* restrict GUy = scratch + nx_all;

  for (int k = 0; k < nfield; ++k) {
    for (int iy = 1; iy < ny_all - 1; ++iy) {
      /* Why does this start at 1 and end at ny-1? Because we need to calculate a
      derivative! Remember, the first and last ng rows/columns of the grid are
      “ghost cells”. The ghost cells exist so that we can implement the periodic
      boundary conditions. We really only care about the canonical cells.

      Suppose that we want to calculate the x derivative of some quantity at the
      (i,j) cell. For this to work, there needs to be (i-1, j) and (i+1, j) cells.
      Thus, i can NOT be 0 or nx_all - 2 (the largest index in cell grid).

      The same basic argument holds for y (we can’t calculate a y derivative if
      y = 0 or y = ny_all - 2). Thus, we only calculate at the cells whose i,j
      indicies are in {1,2.... nx_all-2} x {1,2.. ny_all-2}

      This means that we don’t calculate the derivative on the first/last
      row/column, but this is fine because those rows/columns are ghost cells! */

      int offset = (k*ny_all + iy)*nx_all + 1;
      /* What's offset?
      Remember that F, G, u and v are nfield by nx_all by ny_all arrays.
      For each k in {1,2... n_field}, FU + k*nx_all*ny_all is the starting address
      of the kth subarray of FU. This sub array is stored in ROW MAJOR order.
      Thus, offset is the starting address of the iyth row of the kth sub matrix of
      FU or GU or U or U_half */

      limited_deriv1(FUx + 1, FU + offset, nx_all - 2);
      limited_derivk(GUy + 1, GU + offset, nx_all - 2, nx_all);

      for (int ix = 1; ix < nx_all - 1; ++ix) {
        int offset = (k*ny_all + iy)*nx_all + ix;
        /* What is this offset?
        This is the address of the (ix, iy) entry of the kth subarray of FU/GU/
        U/U_half (remember, there are ghost cells!) */

        // Calculate the (ix, iy) component of the kth component of U_half!
        U_half[offset] = U[offset] - dtcdx2*FUx[ix] - dtcdy2*GUy[ix];
      } // for (int ix = 1; ix < nx_all - 1; ++ix) {
    } // for (int iy = 1; iy < ny_all - 1; ++iy) {
  } // for (int k = 0; k < nfield; ++k) {
} // void central2d_predict(float* restrict U_half,



// Corrector
static
void central2d_correct_sd(float* restrict s,
                          float* restrict d,
                          const float* restrict Uk_x,
                          const float* restrict Uk_y,
                          const float* restrict Uk,
                          const float* restrict FUk,
                          const float* restrict GUk,
                          float dtcdx2,
                          float dtcdy2,
                          int xlo,
                          int xhi)
{
  /* Description: This function calculuates s and d, which are used to
  updated U_half in central2d_correct, for a row of the grid. We calculate s
  and d for each i in {xlo, xlo + 1,... xhi - 1}.

  What are the arguments?
  s, d - values which are used to calculate the "corrected" value of U_half.

  Uk_x - the x derivative of the kth field/component/subarray of U.

  Uk_y - the y derivative of the kth field/component/subarray of U.

  Uk, FUk, GUk - the addresses of a row in the kth field/component/subarrays of
  U, FU, and GK.

  dtcdx2 - (1/2)(dt/dx) (see central2d_step)

  dtcdy2 - (1/2)(dt/dy) (see central2d_step)

  xhi, xlo - the upper and lower indicies within a particular row of the cell
  grid upon which we want to calculate s and d. */

  for (int ix = xlo; ix < xhi; ++ix) {
    s[ix] =
      0.2500f * (Uk [ix] + Uk [ix + 1]) +
      0.0625f * (Uk_x[ix] - Uk_x[ix + 1]) +
      dtcdx2  * (FUk [ix] - FUk [ix + 1]);
  } // for (int ix = xlo; ix < xhi; ++ix) {

  for (int ix = xlo; ix < xhi; ++ix) {
    d[ix] =
      0.0625f * (Uk_y[ix] + Uk_y[ix + 1]) +
      dtcdy2  * (GUk [ix] + GUk [ix + 1]);
  } // for (int ix = xlo; ix < xhi; ++ix) {
} // void central2d_correct_sd(float* restrict s,



// Corrector
static
void central2d_correct(float* restrict U_half,
                       float* restrict scratch,
                       const float* restrict U,
                       const float* restrict FU,
                       const float* restrict GU,
                       float dtcdx2,
                       float dtcdy2,
                       int xlo,
                       int xhi,
                       int ylo,
                       int yhi,
                       int nx_all,
                       int ny_all,
                       int nfield)
{
  /* Description: This function "corrects" U_half (U on a grid which is
  staggered by half a step with respect to the grid that U is on). The values of
  U_half from central2d_predict were used to
  update FU and GU. We then use these values to "correct" or update the values
  of U_half.

  What are the arguments?
  U_half, stratch, U, FU, GU - members (of the same name) of the central2d
  structure.

  dtcdx2 - (1/2)(dt/dx) (see central2d_step).

  dtcdy2 - (1/2)(dt/dy) (see central2d_step).

  xlo - the index of the first column of canonical cells

  xhi - the index of the last columns of canonical cells

  ylo - the index of the first row of canonical cells

  yhi - the index of the last row of canonical cells.

  nx_all - number of columns (x direction) of cells (including ghost cells).

  ny_all - number of rows (y direction) of cells (including ghost cells).

  nfield - the number of fields/subarrays/components in U, U_half, FU, FG. */

  assert(0 <= xlo && xlo < xhi && xhi <= nx_all);
  assert(0 <= ylo && ylo < yhi && yhi <= ny_all);


  // these hold the derivatives of u in the x and y directions.
  float* restrict U_x = scratch;
  float* restrict U_y = scratch +   nx_all;

  float* restrict s0 = scratch + 2*nx_all;
  float* restrict d0 = scratch + 3*nx_all;
  float* restrict s1 = scratch + 4*nx_all;
  float* restrict d1 = scratch + 5*nx_all;

  for (int k = 0; k < nfield; ++k) {

    // Thus, U_half_k/Uk/FUk/GUk are the starting addresses of the kth sub-arrays
    // of U_half/U/FU/GU.
    float*       restrict U_half_k = U_half + k*ny_all*nx_all;
    const float* restrict Uk  = U  + k*ny_all*nx_all;
    const float* restrict FUk = FU + k*ny_all*nx_all;
    const float* restrict GUk = GU + k*ny_all*nx_all;

    // Calculate derivatives of U in the x and y directions, use them to find
    // s and d.
    limited_deriv1(U_x + 1, Uk + ylo*nx_all + 1, nx_all - 2);
    limited_derivk(U_y + 1, Uk + ylo*nx_all + 1, nx_all - 2, nx_all);
    central2d_correct_sd( s1,
                          d1,
                          U_x,                             // Uk_x
                          U_y,                             // Uk_y
                          Uk + ylo*nx_all,                 // Uk
                          FUk + ylo*nx_all,                // FUk
                          GUk + ylo*nx_all,                // GUk
                          dtcdx2,
                          dtcdy2,
                          xlo,
                          xhi);

    for (int iy = ylo; iy < yhi; ++iy) {
      // swap s0 and s1... swap d0 and d1.
      float* tmp;
      tmp = s0; s0 = s1; s1 = tmp;
      tmp = d0; d0 = d1; d1 = tmp;

      // calculate derivatives of U in the x and y directions, use these to
      // calculate s and d.
      limited_deriv1(U_x + 1, Uk + (iy + 1)*nx_all + 1, nx_all - 2);
      limited_derivk(U_y + 1, Uk + (iy + 1)*nx_all + 1, nx_all - 2, nx_all);
      central2d_correct_sd( s1,
                            d1,
                            U_x,                           // Uk_x
                            U_y,                           // Uk_y
                            Uk + (iy + 1)*nx_all,          // Uk
                            FUk + (iy + 1)*nx_all,         // FUk
                            GUk + (iy + 1)*nx_all,         // GUk
                            dtcdx2,
                            dtcdy2,
                            xlo,
                            xhi);

      // Update U_half.
      for (int ix = xlo; ix < xhi; ++ix) {
        U_half_k[ix + iy*nx_all] = (s1[ix] + s0[ix]) - (d1[ix] - d0[ix]);
      } //  for (int ix = xlo; ix < xhi; ++ix) {
    } // for (int iy = ylo; iy < yhi; ++iy) {
  } // for (int k = 0; k < nfield; ++k) {
} // void central2d_correct(float* restrict U_half,



static
void central2d_step(float* restrict U,
                    float* restrict U_half,
                    float* restrict scratch,
                    float* restrict FU,
                    float* restrict GU,
                    int io,
                    int nx,
                    int ny,
                    int ng,
                    int nfield,
                    flux_t flux,
                    speed_t speed,
                    float dt,
                    float dx,
                    float dy)
{
  /* Description: This function completes a time step.

  What are the arguments?
  U_half, U, stratch, FU, GU - the members (of the same name) of the central2d
  structure.

  io - (time step)% 2

  nx - number of canonical cells in the x direction.

  ny - number of canonical cells in the y direction.

  ng - number of layers of ghost cells.

  nfield - number of fields/components/subarrays of FU/GU/U/U_half.

  flux, speed - pointers to the flux and speed functions (see see shallow2d.c)

  dt - time step

  dx, dy - dimensions of a cell. */

  int nx_all = nx + 2*ng;
  int ny_all = ny + 2*ng;

  float dtcdx2 = 0.5*(dt/dx);
  float dtcdy2 = 0.5*(dt/dy);

  flux(FU, GU, U, nx_all*ny_all, nx_all*ny_all);

  // Predictor step.
  central2d_predict(U_half,
                    scratch,
                    U,
                    FU,
                    GU,
                    dtcdx2,
                    dtcdy2,
                    nx_all,
                    ny_all,
                    nfield);

  // Calculate FU and GU using the values of U from the predictor step.
  for (int iy = 1; iy < ny_all - 1; ++iy) {
      int jj = iy*nx_all + 1;
      flux(FU + jj, GU + jj, U_half + jj, nx_all - 2, nx_all*ny_all);
  } // for (int iy = 1; iy < ny_all - 1; ++iy) {

  // Corrector step.
  central2d_correct(U_half + io*(nx_all + 1),
                    scratch,
                    U,
                    FU,
                    GU,
                    dtcdx2,
                    dtcdy2,
                    ng - io,           // xlo
                    nx + ng - io,      // xhi
                    ng - io,           // ylo
                    ny + ng - io,      // yhi
                    nx_all,
                    ny_all,
                    nfield);
} // void central2d_step(float* restrict U,...


/**
 * ### Advance a fixed time
 *
 * The `run` method advances from time 0 (initial conditions) to time
 * `tfinal`.  Note that `run` can be called repeatedly; for example,
 * we might want to advance for a period of time, write out a picture,
 * advance more, and write another picture.  In this sense, `tfinal`
 * should be interpreted as an offset from the time represented by
 * the simulator at the start of the call, rather than as an absolute time.
 *
 * We always take an even number of steps so that the solution
 * at the end lives on the main grid instead of the staggered grid.
 */

static
int central2d_xrun(float* restrict U,
                   float* restrict U_half,
                   float* restrict scratch,
                   float* restrict FU,
                   float* restrict GU,
                   int nx,
                   int ny,
                   int ng,
                   int nfield,
                   flux_t flux,
                   speed_t speed,
                   float tfinal,
                   float dx,
                   float dy,
                   float cfl,
                   float* restrict U_global,
                   float* restrict shared_buffer,
                   int nx_global,
                   int ny_global,
                   const int xlow_local,
                   const int ylow_local)
{
  /* Description: This function runs a simulation!

  What are the arguments?
  U, U_half, FU, GU - members (of the same name) of the central2d structure for
  a piece of the partition of the global array.

  nx - number of columns of canonical cells in U.

  ny  - number of rows of canonical cells in U.

  ng - number of layers of ghost cells.

  nfield - number of fields/components/subarrays of FU/GU/U/U_half. In our case,
  this will be 3.

  flux - pointer to a function which will updated F and G using U (see
  shallow2d.c)

  speed - a function which will find the maximum wave velocity in the x and y
  directions (see shallow2d.c).

  tfinal - the final time we want to work twoards.

  dx, dy - dimensions of a cell.

  cfl - used to determine the time step.

  U_global - a pointer to the U array of the global sim.

  shared_buffer - a pointer to a shared array. This array must have at least
  1 element . This is used fo faciliate communication between the threads. In
  particular, we will use this to calculate dt.

  nx_global - number of rows of canonical cells in U_global.

  ny_global - number of rows of canonical cells in U_global.

  xlow_local - the rows in sim_local's U array correspond to rows xlow_local to
  xlow_local + nx in U.

  ylow_local - the columns in sim_local's U array correspond to columns
  ylow_local to ylow_local + ny in U.*/

  // Set up for the main loop!
  int nstep = 0;
  int nx_all = nx + 2*ng;
  int ny_all = ny + 2*ng;
  bool done = false;
  float t = 0;

  while (!done) {
    float cxy[2] = {1.0e-15f, 1.0e-15f};
    /* What's going on here ( ^ )?
    We set the elements of cxy to small but non-zero values so that when we
    calculate dt, we don’t divide by zero! */

    /* Apply periodic boundary conditions to U_global. This ensures that the
    ghost cells of U_global can easily be pulled into the ghost cells of each
    processors U. */

    /* one way to improve this would be to eliminate the single directive.
    Throw in for directives inside of this function, so that the different
    threads complete different parts of the BCs in parallel. */
    #pragma omp single
    {
      central2d_periodic( U_global,
                          nx_global,
                          ny_global,
                          ng,
                          nfield);
    } // #pragma omp single

    /* Get boundary information for U from U_global */
    central2d_local_BC(U,
                       nx,
                       ny,
                       ng,
                       nfield,
                       U_global,
                       nx_global,
                       ny_global,
                       xlow_local,
                       ylow_local);

    /* Calculate maximum wave speed in the x and y directions, use this and cfl
    to determine dt.

    To do this, each thread calculates the maximum x and y velocity for
    it's U and uses this calculate dt given its data. The threads then one-by-
    one store their values into the shared buffer.

    Once every thread has written its data to the shared buffer, the buffer
    contains the minimum dt, which each processor can then use. This is
    essentially doing what an MPI reduce operation would do, but using openMP
    to do it. */

    /* calculate dt using the data on our piece of the partition. */
    speed(cxy, U, nx_all*ny_all, nx_all*ny_all);
    float dt_local = cfl / fmaxf(cxy[0]/dx, cxy[1]/dy);

    /* The first processor to complete its work can initialize the
    shared buffer (we need to do this, the buffer may contain 0, which would
    wreak havoc on the rest of the code). */
    #pragma omp single
    {
      shared_buffer[0] = dt_local;
    } // #pragma omp single

    // Now, each thread writes its data to shared_buffer one by one.
    #pragma omp critical
    {
      shared_buffer[0] = fmin(dt_local, shared_buffer[0]);
    } // #pragma omp critical
    #pragma omp barrier

    /* shared_buffer[0] now contains the minimum dt (across the threads), we
    can read it in to move on. */
    float dt = shared_buffer[0];

    /* Check if we are ready to stop looping. This is how the loop eventually
    stops and ensures that we always stop at t final (think about it). */
    if (t + 2*dt >= tfinal) {
      dt = (tfinal - t)/2;
      done = true;
    } // if (t + 2*dt >= tfinal) {

    // Compute an odd time step.
    central2d_step( U,
                    U_half,
                    scratch,
                    FU,
                    GU,
                    0,                 // io
                    nx + 4,            // nx
                    ny + 4,            // ny
                    ng - 2,            // ng
                    nfield,
                    flux,
                    speed,
                    dt,
                    dx,
                    dy);

    // Compute an even time step. Note that the roles of U and U_half
    // are swapped.
    central2d_step( U_half,
                    U,
                    scratch,
                    FU,
                    GU,
                    1,                 // io
                    nx,
                    ny,
                    ng,
                    nfield,
                    flux,
                    speed,
                    dt,
                    dx,
                    dy);
    t += 2*dt;
    nstep += 2;

    // copy boundary of U to corresponding entries of U global
    central2d_local_to_global(U,
                              nx,
                              ny,
                              ng,
                              nfield,
                              U_global,
                              nx_global,
                              ny_global,
                              xlow_local,
                              ylow_local);

    /* Wait until all threads are done copying their data into U. We need
    to do this because at the start of the next time step, we're going to
    apply periodic BCs, and that will only give the result we want if each
    thread has written its cells into U_global. */
    #pragma omp barrier
  } // while (!done) {

  // return the number of time steps.
  return nstep;
} // int central2d_xrun(float* restrict U,



int central2d_run(central2d_t* sim_local,
                  central2d_t* sim,
                  const int xlow_local,
                  const int ylow_local,
                  float tfinal)
{
  /* Description: This is a wrapper for central2d_xrun.

  What are the arguments?
  sim_local - a central2d structure which corresponds to a piece of the global
  grid.

  sim - a central2d strcuture which corresponds to the global grid.

  U_global - a pointer to the U array of the global simulation

  xlow_local - the rows in sim_local's U array correspond to rows xlow_local to
  xlow_local + sim_local->nx in U.

  ylow_local - the columns in sim_local's U array correspond to columns
  ylow_local to ylow_local + sim_local->ny in U.

  tfinal - the final time in the simulation */

  return central2d_xrun(sim_local->U,
                        sim_local->U_half,
                        sim_local->scratch,
                        sim_local->FU,
                        sim_local->GU,
                        sim_local->nx,
                        sim_local->ny,
                        sim_local->ng,
                        sim_local->nfield,
                        sim_local->flux,
                        sim_local->speed,
                        tfinal,
                        sim_local->dx,
                        sim_local->dy,
                        sim_local->cfl,
                        sim->U,                  // U_global
                        sim->scratch,            // shared_buffer
                        sim->nx,                 // nx_global
                        sim->ny,                 // ny_global
                        xlow_local,
                        ylow_local);
} // int central2d_run(central2d_t* sim_local,...

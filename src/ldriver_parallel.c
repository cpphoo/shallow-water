#include "stepper_parallel.h"
#include "shallow2d.h"

#ifdef _OPENMP
#include <omp.h>
#elif defined SYSTIME
#include <sys/time.h>
#endif

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <assert.h>
#include <stdio.h>

//ldoc on
/**
 * # Driver code
 *
 * The driver code is where we put together the time stepper and
 * the physics routines to actually solve the equations and make
 * pretty pictures of the solutions.
 *
 * ## Diagnostics
 *
 * The numerical method is supposed to preserve (up to rounding
 * errors) the total volume of water in the domain and the total
 * momentum.  Ideally, we should also not see negative water heights,
 * since that will cause the system of equations to blow up.  For
 * debugging convenience, we'll plan to periodically print diagnostic
 * information about these conserved quantities (and about the range
 * of water heights).
 */

void solution_check(central2d_t* sim)
{
  int nx = sim->nx, ny = sim->ny;
  float* U = sim->U;
  float h_sum = 0, hu_sum = 0, hv_sum = 0;
  float hmin = U[central2d_offset(sim,0,0,0)];
  float hmax = hmin;

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      float h = U[central2d_offset(sim,0,i,j)];
      h_sum += h;
      hu_sum += U[central2d_offset(sim,1,i,j)];
      hv_sum += U[central2d_offset(sim,2,i,j)];
      hmax = fmaxf(h, hmax);
      hmin = fminf(h, hmin);
    } // for (int i = 0; i < nx; ++i) {
  } // for (int j = 0; j < ny; ++j) {

  float cell_area = sim->dx * sim->dy;
  h_sum *= cell_area;
  hu_sum *= cell_area;
  hv_sum *= cell_area;
  printf("-\n  Volume: %g\n  Momentum: (%g, %g)\n  Range: [%g, %g]\n",
         h_sum, hu_sum, hv_sum, hmin, hmax);
  assert(hmin > 0);
} // void solution_check(central2d_t* sim)

/**
 * ## I/O
 *
 * After finishing a run (or every several steps), we might want to
 * write out a data file for further processing by some other program
 * -- in this case, a Python visualizer.  The visualizer takes the
 * number of pixels in x and y in the first two entries, then raw
 * single-precision raster pictures.
 */

FILE* viz_open(const char* fname, central2d_t* sim, int vskip)
{
  FILE* fp = fopen(fname, "w");

  if (fp) {
    float xy[2] = {sim->nx/vskip, sim->ny/vskip};
    fwrite(xy, sizeof(float), 2, fp);
  } // if (fp) {

  return fp;
} // FILE* viz_open(const char* fname, central2d_t* sim, int vskip)

void viz_close(FILE* fp)
{
  fclose(fp);
} // void viz_close(FILE* fp)

void viz_frame(FILE* fp, central2d_t* sim, int vskip)
{
  if (!fp) {
    return;
  } // if (!fp) {

  for (int iy = 0; iy < sim->ny; iy += vskip) {
    for (int ix = 0; ix < sim->nx; ix += vskip) {
      fwrite(sim->U + central2d_offset(sim,0,ix,iy),
             sizeof(float), 1, fp);
    } // for (int ix = 0; ix < sim->nx; ix += vskip) {
  } // for (int iy = 0; iy < sim->ny; iy += vskip) {
} // void viz_frame(FILE* fp, central2d_t* sim, int vskip)

/**
 * ## Lua driver routines
 *
 * A better way to manage simulation parameters is by a scripting
 * language.  Python is a popular choice, but I prefer Lua for many
 * things (not least because it is an easy build).  It's also quite
 * cheap to call a Lua function for every point in a mesh
 * (less so for Python, though it probably won't make much difference).
 *
 * ### Lua callback functions
 *
 * We specify the initial conditions by providing the simulator
 * with a callback function to be called at each cell center.
 * The callback function is assumed to be the `init` field of
 * a table at index 1.
 */

void lua_init_sim(lua_State* L, central2d_t* sim)
{
  lua_getfield(L, 1, "init");
  if (lua_type(L, -1) != LUA_TFUNCTION) {
    luaL_error(L, "Expected init to be a string");
  } // if (lua_type(L, -1) != LUA_TFUNCTION) {

  int nx = sim->nx, ny = sim->ny, nfield = sim->nfield;
  float dx = sim->dx, dy = sim->dy;
  float* U = sim->U;

  for (int ix = 0; ix < nx; ++ix) {
    float x = (ix + 0.5) * dx;

    for (int iy = 0; iy < ny; ++iy) {
      float y = (iy + 0.5) * dy;
      lua_pushvalue(L, -1);
      lua_pushnumber(L, x);
      lua_pushnumber(L, y);
      lua_call(L, 2, nfield);

      for (int k = 0; k < nfield; ++k) {
        U[central2d_offset(sim,k,ix,iy)] = lua_tonumber(L, k-nfield);
      } // for (int k = 0; k < nfield; ++k) {

      lua_pop(L, nfield);
    } // for (int iy = 0; iy < ny; ++iy) {
  } // for (int ix = 0; ix < nx; ++ix) {

  lua_pop(L,1);
} // void lua_init_sim(lua_State* L, central2d_t* sim)


/**
 * ### Running the simulation
 *
 * The `run_sim` function looks a lot like the main routine of the
 * "ordinary" command line driver.  We specify the initial conditions
 * by providing the simulator with a callback function to be called at
 * each cell center.  Note that we have two different options for
 * timing the steps -- we can use the OpenMP timing routines
 * (preferable if OpenMP is available) or the POSIX `gettimeofday`
 * if the `SYSTIME` macro is defined.  If there's no OpenMP and
 * `SYSTIME` is undefined, we fall back to just printing the number
 * of steps without timing information.
 */

int run_sim(lua_State* L)
{
  // lua set up.
  int n = lua_gettop(L);
  if (n != 1 || !lua_istable(L, 1)) {
    luaL_error(L, "Argument must be a table");
  } // if (n != 1 || !lua_istable(L, 1)) {

  lua_getfield(L, 1, "w");
  lua_getfield(L, 1, "h");
  lua_getfield(L, 1, "cfl");
  lua_getfield(L, 1, "ftime");
  lua_getfield(L, 1, "nx");
  lua_getfield(L, 1, "ny");
  lua_getfield(L, 1, "vskip");
  lua_getfield(L, 1, "frames");
  lua_getfield(L, 1, "out");

  double grid_width     = luaL_optnumber( L, 2,  2.0);
  double grid_height    = luaL_optnumber( L, 3,  grid_width);
  double cfl            = luaL_optnumber( L, 4,  0.45);
  double ftime          = luaL_optnumber( L, 5,  0.01);
  int nx                = luaL_optinteger(L, 6,  200);
  int ny                = luaL_optinteger(L, 7,  nx);
  int vskip             = luaL_optinteger(L, 8,  1);
  int frames            = luaL_optinteger(L, 9,  50);
  const char* fname     = luaL_optstring( L, 10, "sim.out");
  lua_pop(L, 9);

  // initialize the global simulation object (which holds the entire grid)
  central2d_t* sim = central2d_init(grid_width,
                                    grid_height,
                                    nx,
                                    ny,
                                    3,                     // nfield
                                    shallow2d_flux,
                                    shallow2d_speed,
                                    cfl);

  // more lua stuff.
  lua_init_sim(L,sim);
  printf("%g %g %d %d %g %d %g\n", grid_width, grid_height, nx, ny, cfl, frames, ftime);
  FILE* viz = viz_open(fname, sim, vskip);
  solution_check(sim);
  viz_frame(viz, sim, vskip);


  //////////////////////////////////////////////////////////////////////////////
  // Begin Parallel region!
  const int n_rows = 1;
  const int n_cols = 2;

  double tcompute = 0;

  // First, make sure that openMP is defined... if not, then abort!
  #ifndef _OPENMP
    printf("openMP not defined. Aborting\n");
    abort();
  #endif

  #pragma omp parallel num_threads(n_rows*n_cols)
  {
    // First, check that there are n_rows*n_cols threads
    if(omp_get_num_threads() != n_rows*n_cols) {
      printf("Couldn't create enough threads! Wanted %d but got %d\n aborting.\n", n_rows*n_cols, omp_get_num_threads());
      abort();
    } // if(omp_get_num_threads() != n_rows*n_cols) {

    /* In this approach, we partition the global grid into a bunch of
    sub grids. We split the rows of the global grid into n_rows rows and
    n_cols cols. We assign threads to pieces of the partition based on their
    thread number.

    Let's do that now. pieces of the partition are indentified by two indicies,
    p_row and p_col. These specify the row and column (within the partition) of the
    piece. */
    int p_row = omp_get_thread_num() % n_rows;
    int p_col = ((int) omp_get_thread_num()) / ((int) n_rows);

    /* Now that each processor has its sub grid indicies, we can determine which
    rows and columns (within the global grid) each processor will be responsible
    for. */
    const int ny_p = ny / n_rows;
    int ylow_local = p_row*ny_p;
    int yhigh_local;
    if(p_row == n_rows) { yhigh_local = ny; }
    else { yhigh_local = ylow_local + ny_p; }

    const int nx_p = nx / n_cols;
    int xlow_local = p_col*nx_p;
    int xhigh_local;
    if(p_col == n_cols) { xhigh_local = nx; }
    else { xhigh_local = xlow_local + nx_p; }

    // Now set up the local simulation structure.
    const int nx_local = xhigh_local - xlow_local;
    const int ny_local = yhigh_local - ylow_local;
    central2d_t* sim_local = central2d_init( grid_width,  // This is wrong, but it's okay... see below.
                                             grid_height, // This is wrong, but it's okay... see below.
                                             nx_local,
                                             ny_local,
                                             3,           // nfield
                                             shallow2d_flux,
                                             shallow2d_speed,
                                             cfl);

    /* Why do we pass the global grid_width and grid_height to the initializer
    for sim_local? Remember, sim_local only deals with a piece of the global
    grid, so its height and width will be smaller. If we look at the code for
    the initializer, we'll notice that grid_width and grid_height are only used
    to calculate dx and dy. Every cell in the global grid has the same size.
    Thus, dx and dy for each local grid should be equal to that of the global
    one. When we initialized the global sim variable, it calculated dx and dy.
    Thus, dx and dy are already known, we just need to get them to sim_local.

    The idea here is to just pass some junk values to central2d_init when
    initializing sim_local. This initializer will calculate incorrect values for
    dx and dy (for the local grid). After initialization is done, we will
    overwrite the faulty values of dx and dy using the ones in sim. */
    sim_local->dx = sim->dx;
    sim_local->dy = sim->dy;

    /* Now that sim_local has been initialized, we need to set up its U array.
    To do that, we need to copy the corresponding part of the global U array
    into the local U array. */
    float* U_local = sim_local->U;
    float* U = sim -> U;
    for(int k = 0; k < 3; ++k) {  // 3 = nfield
      for(int iy = 0; iy < ny_local; ++iy) {
        for(int ix = 0; ix < nx_local; ++ix) {
          /* We need to copy a piece of the global U to the local U.
          The first row of U_local corresponds to row xlow_local in U.

          Likewise, the first column of U_local corresponds to row ylow_local
          in U. */
          U_local[central2d_offset(sim_local, k, ix, iy)] = U[central2d_offset(sim, k, xlow_local + ix , ylow_local + iy)];
        } // for(int ix = 0; ix < sim->nx; ++ix) {
      } // for(int iy = 0; iy < sim->ny; ++iy) {
    } // for(int k = 0; k < 3; ++k) {

    // wait for all threads to set up their local arrays.
    #pragma omp barrier

    /* Now, at long last, the local simulation structures are set up and ready
    to go! Let's cycle through the timesteps! */
    for (int i = 0; i < frames; ++i) {
      double t0 = omp_get_wtime();
      int nstep = central2d_run(sim_local,
                                sim,
                                xlow_local,
                                ylow_local,
                                ftime);
      double t1 = omp_get_wtime();
      double elapsed = t1 - t0;

      /* IO: At the end of run, each thread writes its local U to global U.
      Before we can run IO, we need each processor to be done writing to global.
      Thus, we have a barrier.

      Once everything is synchronized (every thread has written to global U),
      we can print out diagostics on sim and write a frame of sim to file. */
      #pragma omp barrier
      #pragma omp sections
      {
        // One section to run out diagnostic on U.
        #pragma omp section
        {
          solution_check(sim);
          tcompute += elapsed;
          printf("  Time: %e (%e for %d steps)\n", elapsed, elapsed/nstep, nstep);
        } // #pragma omp section

        // One section to write a frame of U to memory.
        #pragma omp section
        {
          viz_frame(viz, sim, vskip);
        } // #pragma omp section
      } // #pragma omp single {
    } // for (int i = 0; i < frames; ++i) {

    // Free the local sim structure.
    central2d_free(sim_local);
  } // #pragma omp parallel num_threads(4)

  printf("Total compute time: %e\n", tcompute);
  viz_close(viz);
  central2d_free(sim);
  return 0;
} // int run_sim(lua_State* L)


/**
 * ### Main
 *
 * The main routine has the usage pattern
 *
 *     lshallow tests.lua args
 *
 * where `tests.lua` has a call to the `simulate` function to run
 * the simulation.  The arguments after the Lua file name are passed
 * into the Lua script via a global array called `args`.
 */

int main(int argc, char** argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s fname args\n", argv[0]);
    return -1;
  } // if (argc < 2) {

  lua_State* L = luaL_newstate();
  luaL_openlibs(L);
  lua_register(L, "simulate", run_sim);

  lua_newtable(L);
  for (int i = 2; i < argc; ++i) {
    lua_pushstring(L, argv[i]);
    lua_rawseti(L, 1, i-1);
  } // for (int i = 2; i < argc; ++i) {
  lua_setglobal(L, "args");

  if (luaL_dofile(L, argv[1])) {
    printf("%s\n", lua_tostring(L,-1));
  } // if (luaL_dofile(L, argv[1])) {
  lua_close(L);
  return 0;
} // int main(int argc, char** argv)

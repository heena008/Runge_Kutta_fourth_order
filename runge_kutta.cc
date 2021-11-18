// Deal.ii
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

// STL
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>

#include <deal.II/base/logstream.h>

#include <deal.II/base/time_stepping.h>

#include "advection_field.hpp"

using namespace dealii;

template <int dim> class RK4 : public Tensor<1, dim> {
public:
  RK4() : Tensor<1, dim>() {}
  virtual Tensor<1, dim> value(const Point<dim> &p, double time,
                               double time_step) const;

private:
  AdvectionField<dim> advection_field;
  double T_max;
};

template <int dim>
Tensor<1, dim> RK4<dim>::value(const Point<dim> &p, double time,
                               double time_step) const {

  Tensor<1, dim> yn;

  /* Here the equation is dy/dt =c; where c is velocity y is the vertex of mesh
   * and t is the time
   *
   */

  auto k1 = time_step * advection_field.value(p),
       k2 = time_step * (advection_field.value(p) + k1 / 2),
       k3 = time_step * (advection_field.value(p) + k2 / 2),
       k4 = time_step * (advection_field.value(p) + k3);
  return yn = (1 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
}
///////////////////////////////////////////////////////////////////////////////////////////////

template <int dim> class MeshDeformer {
public:
  MeshDeformer();

  void run();

  double diffOfy(const double time, const Vector<double> time_new) {
    return time; // function x^2 + y^2
  }

private:
  void make_grid();

  MPI_Comm mpi_communicator;

  /*!
   * Distributed triangulation
   */
  parallel::distributed::Triangulation<dim> triangulation;

  GridOut grid_out;
  AdvectionField<dim> advection_field;
  Vector<double> solution;
  RK4<dim> rk4;

  double time;
  double time_step;
  /*!
   * Final simulation time.
   */
  double T_max;
};

template <int dim>
MeshDeformer<dim>::MeshDeformer()
    : mpi_communicator(MPI_COMM_WORLD),
      triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                        Triangulation<dim>::smoothing_on_refinement |
                        Triangulation<dim>::smoothing_on_coarsening)),
      time(0.0), time_step(0.1), T_max(1.0) {}

template <int dim> void MeshDeformer<dim>::make_grid() {
  GridGenerator::hyper_cube(triangulation, 0, 1, /* colorize faces */
                            false);
  triangulation.refine_global(1);

  while (T_max >= time)
  {

    for (const auto &cell : triangulation.active_cell_iterators()) {

      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; i++) {
        Point<2> &v = cell->vertex(i); // initial vertex

        auto v0 = v;
        // advection_field.set_time(time+time_step);
        char method = 'A';

        switch (method) {
        case 'A':
          v += rk4.value(v, time, time_step); // Runge_kutta_fourth_order
          break;

        case 'B':
          v += advection_field.value(v) * time_step; // foward Euler
          break;
        default:
          printf("Invalid method\n");
        }

        std::cout << v << "  t  " << time << "  dt " << T_max << std::endl;
        std::cout << v0 << std::endl;
        std::cout << advection_field.value(v) << std::endl;
        std::cout << "----------------------------------" << std::endl;

        GridOut grid_out;
        std::ofstream output("mesh-" + std::to_string(time) + ".vtu");
        grid_out.write_vtu(triangulation, output);

        std::ofstream out("grid" + std::to_string(time) + ".svg");
        grid_out.write_svg(triangulation, out);
        std::cout << "Grid written to grid-2.svg" << std::endl;
        time = time + time_step;
      }
    }
  }
}

template <int dim> void MeshDeformer<dim>::run()
{
	make_grid();
}

int main(int argc, char *argv[])
{
  try {
    Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
    deallog.depth_console(0);
    MeshDeformer<2> rk4_problem_2d;
    rk4_problem_2d.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  return 0;
}

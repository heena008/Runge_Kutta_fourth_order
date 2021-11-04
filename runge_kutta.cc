// Deal.ii
#include <deal.II/base/tensor_function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

// STL
#include <iostream>
#include <fstream>
#include <cmath>


//#include "reconstruction_mesh.hpp"

#include <deal.II/base/logstream.h>

using namespace dealii;


template <int dim> class AdvectionField : public TensorFunction<1, dim> {
public:
  AdvectionField() : TensorFunction<1, dim>() {}

  virtual Tensor<1, dim> value(const Point<dim> &point) const override;
  virtual void value_list(const std::vector<Point<dim>> &points,
                          std::vector<Tensor<1, dim>> &values) const override;

private:
  const double pi = numbers::PI;
};

template <int dim>
Tensor<1, dim> AdvectionField<dim>::value(const Point<dim> &p) const
{
  const double t = this->get_time();

  Tensor<1, dim> value;
  value.clear();

  double t0 = 0;

    const double tn = 1;

    const double n=10;

    const double dt = (tn-t0)/n;

    const double c=3;

  // Here velocity is consider to be constant let it be 1 m/s.

  for (unsigned int d = 0; d < dim; ++d)

    {


      value[d] =0.8 * std::sin( 3.14* p[1]);

     // value[d] =1; //0.8 * std::sin( numbers::PI * p[1]);

    }

  return value;
}

template <int dim>
void AdvectionField<dim>::value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Tensor<1, dim>> &values) const {
  Assert(points.size() == values.size(),
         ExcDimensionMismatch(points.size(), values.size()));

  for (unsigned int p = 0; p < points.size(); ++p) {
    values[p].clear();
    values[p] = value(points[p]);
  }
}

///////////////////////////////////////////////////////////////

template <int dim> class RK4 : public Tensor<1, dim> {
public:
  RK4() : Tensor<1, dim>() {}
  virtual Tensor<1, dim> value(const Point<dim> &p) const;

private:
  AdvectionField<dim> advection_field;
};

template <int dim> Tensor<1, dim> RK4<dim>::value(const Point<dim> &p) const {

  Tensor<1, dim> yn;

  double t0 = 0;

  const double tn = 1;

  const double n=10;

  const double dt = (tn-t0)/n;

  /* Here the equation is dy/dt =c; where c is velocity y is the vertex of mesh and t is the time
   *
   */

  yn = advection_field.value(p);

  for (double t0 = 0; t0 <= tn; t0++) {
    auto k1 = dt * advection_field.value(p),
         k2 = dt * (advection_field.value(p) + k1 / 2),
         k3 = dt * (advection_field.value(p) + k2 / 2),
         k4 = dt * (advection_field.value(p) + k3);
    return yn += (1 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    t0 = t0 + dt;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////

template <int dim> class MeshDeformer {
public:
  MeshDeformer();

  void run();

private:
  void make_grid();

  Triangulation<dim> triangulation;

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
    : triangulation(), time(1.0), time_step(1), T_max(1.0) {}

template <int dim> void MeshDeformer<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, 0, 1, /* colorize faces */
                            false);
  triangulation.refine_global(4);

  std::vector<Tensor<1, dim>> vetices_on_cell;
   		 while (T_max>=0.09) {
  for (const auto &cell : triangulation.active_cell_iterators()) {

    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; i++) {
      Point<2> &v = cell->vertex(i); // initial vertex

      rk4.value(v);

      std::cout << "----------------------------------" << std::endl;

      std::cout << rk4.value(v) << "  rk4    " << std::endl;
      std::cout << "----------------------------------" << std::endl;

      GridOut       grid_out;
       											std::ofstream output("mesh-" + std::to_string(T_max) + ".vtu");
       											grid_out.write_vtu( triangulation, output);
       										 T_max=T_max-time_step;
    }
  }
  }

}

template <int dim> void MeshDeformer<dim>::run()
{
	make_grid();
}

int main()
{

  MeshDeformer<2> rk4_problem_2d;
  rk4_problem_2d.run();

  return 0;
}

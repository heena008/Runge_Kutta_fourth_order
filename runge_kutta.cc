// Deal.ii
#include <deal.II/base/tensor_function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

// STL
#include <iostream>
#include <fstream>
#include <cmath>

#include <deal.II/base/logstream.h>

#include "advection_field.hpp"

using namespace dealii;


template <int dim> class CharacteristicEquation : public TensorFunction<1, dim> {
public:
  CharacteristicEquation() : TensorFunction<1, dim>() {}

  virtual Tensor<1, dim> value(const Point<dim> &point) const override;
  virtual void value_list(const std::vector<Point<dim>> &points,
                          std::vector<Tensor<1, dim>> &values) const override;

private:
  const double pi = numbers::PI;

  AdvectionField<dim> advection_field;
};

template <int dim>
Tensor<1, dim> CharacteristicEquation<dim>::value(const Point<dim> &p) const
{
  const double t = this->get_time();

  Tensor<1, dim> value;
  value.clear();

  double t0 = 0;

    const double tn = 1;

    const double n=10;

    const double dt = (tn-t0)/n;



  // Here velocity is consider to be constant let it be 1 m/s.

    const double T = 4;

    for (unsigned int d = 0; d < dim; ++d)

     {
    	auto c=advection_field.value(p);

       value[d] =p[d]+c[d]*dt;

      // value[d] =1; //0.8 * std::sin( numbers::PI * p[1]);

     }
      return value;




  return value;
}


template <int dim>
void CharacteristicEquation<dim>::value_list(
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

template <int dim>
class RK4 : public Tensor<1, dim>
{
public:
  RK4() : Tensor<1, dim>() {}
  virtual Tensor<1, dim> value(const Point<dim> &p, double time, double time_step) const;

private:
  CharacteristicEquation<dim> charac_field;
  double T_max;
};


template <int dim> Tensor<1, dim> RK4<dim>::value(const Point<dim> &p,double time, double time_step) const {

  Tensor<1, dim> yn;

double T_max=1.0;

  /* Here the equation is dy/dt =c; where c is velocity y is the vertex of mesh and t is the time
   *
   */

    auto k1 = time_step * charac_field.value(p),
         k2 = time_step * (charac_field.value(p) + k1 / 2),
         k3 = time_step * (charac_field.value(p) + k2 / 2),
         k4 = time_step * (charac_field.value(p) + k3);
    return yn = (1 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

}
///////////////////////////////////////////////////////////////////////////////////////////////

template <int dim> class MeshDeformer {
public:
  MeshDeformer();

  void run();

private:
  void make_grid();

  MPI_Comm mpi_communicator;

    /*!
     * Distributed triangulation
     */
    parallel::distributed::Triangulation<dim> triangulation;

  GridOut grid_out;
  CharacteristicEquation<dim> charac_field;
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
		                      Triangulation<dim>::smoothing_on_coarsening)), time(0.0), time_step(0.1), T_max(1.0) {}

template <int dim> void MeshDeformer<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, 0, 1, /* colorize faces */
                            false);
  triangulation.refine_global(1);

  std::vector<Tensor<1, dim>> vetices_on_cell;

  while (T_max >= time) {

      for (const auto &cell : triangulation.active_cell_iterators()) {

        for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; i++)
        {
          Point<2> &v = cell->vertex(i); // initial vertex

          auto v0 = v;
          //    							advection.set_time(time-time_step);
          char method = 'A';

             switch(method) {
                case 'A' :
                	v+=rk4.value(v,time,time_step);//Runge_kutta_fourth_order
                   break;

                case 'B' :
                	  v+=charac_field.value(v)*time_step;//foward Euler
                   break;
                default :
                   printf("Invalid method\n" );
             }
         //

//          auto  k1 = time_step * charac_field.value(v),
//           											 		k2 = time_step * ( charac_field.value(v)+ k1 / 2),
//           											 		k3 = time_step * (charac_field.value(v) + k2 / 2),
//           													k4 = time_step * (charac_field.value(v) + k3);
//           											 	v+=(1/6.0)*(k1+2*k2+2*k3+k4 ) ;
          //
          std::cout << v << "  t  "<<time <<"  dt "<< T_max<<std::endl;
//          std::cout << v0 << std::endl;
          std::cout << charac_field.value(v) << std::endl;
          std::cout << "----------------------------------" << std::endl;

          GridOut grid_out;
          std::ofstream output("mesh-" + std::to_string(time) + ".vtu");
          grid_out.write_vtu(triangulation, output);


          std::ofstream out("grid"+ std::to_string(time)+ ".svg" );
             		  			    				  grid_out.write_svg(triangulation, out);
             		  			    				  std::cout << "Grid written to grid-2.svg" << std::endl;

        }

      }
      time =time+time_step;
//      std::ofstream out("grid-2.svg");
//      		  			    				  GridOut       grid_out;
//      		  			    				  grid_out.write_svg(triangulation, out);
//      		  			    				  std::cout << "Grid written to grid-2.svg" << std::endl;
    }


}

template <int dim> void MeshDeformer<dim>::run()
{
	make_grid();
}

int main(int argc, char *argv[])
{
  try
    {
	  Utilities::MPI::MPI_InitFinalize mpi_initialization(
	  	      argc, argv, numbers::invalid_unsigned_int);
	  		deallog.depth_console (0);
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


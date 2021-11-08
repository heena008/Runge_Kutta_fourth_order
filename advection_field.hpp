/*
 * advection_field.hpp
 *
 *  Created on: Nov 4, 2019
 *      Author: heena
 */

#ifndef INCLUDE_ADVECTION_FIELD_HPP_
#define INCLUDE_ADVECTION_FIELD_HPP_

// Deal.ii
#include <deal.II/base/tensor_function.h>

// STL
#include <cmath>
#include <fstream>

using namespace dealii;

/*!
 * @class AdvectionField
 * @brief Class implements Advection Field.
 */

template <int dim> class AdvectionField : public TensorFunction<1, dim>
{
public:
  AdvectionField() : TensorFunction<1, dim>() {}

  virtual Tensor<1, dim> value(const Point<dim> &point) const override;
  virtual void value_list(const std::vector<Point<dim>> &points,
						  std::vector<Tensor<1, dim>> &values) const override;
};

template <int dim>
Tensor<1, dim> AdvectionField<dim>::value(const Point<dim> &p) const
{
  const double t = this->get_time();

  Tensor<1, dim> value;
  value.clear();

  for (unsigned int d = 0; d < dim; ++d)

  {


    value[d] =1;

   // value[d] =1; //0.8 * std::sin( numbers::PI * p[1]);

  }
   return value;

}

template <int dim>
void AdvectionField<dim>::value_list(const std::vector<Point<dim>> &points,
									 std::vector<Tensor<1, dim>> &values) const
{
	Assert(points.size() == values.size(),
	ExcDimensionMismatch(points.size(), values.size()));

  for (unsigned int p = 0; p < points.size(); ++p)
  {

    values[p].clear();
    values[p] = value(points[p]);

  }
}


#endif /* INCLUDE_ADVECTION_FIELD_HPP_ */

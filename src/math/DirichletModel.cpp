/*
 * DirichletModel.cpp
 *
 *  Created on: Jun 16, 2016
 *      Author: zhengqi
 */

#include "DirichletModel.h"

namespace EGriceLab {
namespace Math {
using namespace std;
using namespace Eigen;

const double DirichletModel::DEFAULT_ETA = 0.001;
const double DirichletModel::DEFAULT_ABS_EPS_COST = 1e-6;
const double DirichletModel::DEFAULT_ABS_EPS_PARAMS = 1e-6;
const double DirichletModel::DEFAULT_REL_EPS_COST = 0;
const double DirichletModel::DEFAULT_REL_EPS_PARAMS = 0;
const IOFormat DirichletModel::FULL_FORMAT(Eigen::FullPrecision);

} /* namespace Math */
} /* namespace EGriceLab */



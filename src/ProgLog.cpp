/*
 * ProgLog.cpp
 *
 *  Created on: Dec 14, 2016
 *      Author: zhengqi
 */

#include "ProgLog.h"

namespace EGriceLab {

using namespace std;

ProgLog errorLog(cerr, ProgLog::LOG_ERROR);
ProgLog warningLog(cerr, ProgLog::LOG_WARNING);
ProgLog infoLog(cerr, ProgLog::LOG_INFO);

} /* namespace EGriceLab */

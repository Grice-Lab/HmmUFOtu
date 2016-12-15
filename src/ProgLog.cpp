/*
 * ProgLog.cpp
 *
 *  Created on: Dec 14, 2016
 *      Author: zhengqi
 */

#include "ProgLog.h"

namespace EGriceLab {

using namespace std;

ProgLog errorLog(cerr, LOG_ERROR);
ProgLog warningLog(cerr, LOG_WARNING);
ProgLog infoLog(cerr, LOG_INFO);
ProgLog debugLog(cerr, LOG_DEBUG);

} /* namespace EGriceLab */

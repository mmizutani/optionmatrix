/*
   optionmatrix:

   Options & Futures Matrix Modeler
   View and Control Theoretical Option Chains

   File: futures.c of optionmatrix

   Copyright (c) Anthony Bradford. 2012.
   http://opensourcefinancialmodels.com
   info@opensourcefinancialmodels.com

   optionmatrix may be freely redistributed.
   See file COPYING included with this distribution for license information
*/

/* 
   optionmatrix program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   optionmatrix program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>

#include "basicmodels.h"

double futures_price(const double spot,const double rate, const double dividend, const double time)
{
    if(time<=0)
    {
        return spot;
    }
    /* dividend is the continuous dividend yield */
    
    return exp((rate-dividend)*time)*spot;
}

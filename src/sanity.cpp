/* optionmatrix:                                                            */
/*                                                                          */
/* Options & Futures Matrix Modeler                                         */
/* View and Control Theoretical Option Chains                               */
/*                                                                          */
/* File: sanity.cpp of optionmatrix                                         */
/*                                                                          */
/* Copyright (c) Anthony Bradford. 2012.                                    */
/* http://opensourcefinancialmodels.com                                     */
/* info@opensourcefinancialmodels.com                                       */
/*                                                                          */
/* optionmatrix may be freely redistributed.                                */
/* See file COPYING included with this distribution for license information */

/* 
   optionmatrix is free software; you can redistribute it and/or modify
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

#include "defs.h"
#include "extern.h"
#include "gtk_extern.h"

int sanity_check(struct _properties *properties, char *statusMessage)
{
    bool messageSet = 0;

    if( day_changed() && properties->decimalorcalendar == CALENDAR )
    {
      expires(LEG1,30,properties->skipmonth,properties);
      expires(LEG2,30,properties->skipmonth2,properties);
    }

    if( option_algorithms[properties->modeltype].iUseZ && properties->data.UseZ == 0 &&
        option_algorithms[properties->modeltype].bZallow0Negative == 0 &&
        option_algorithms[properties->modeltype].Zdefault != 0 )
    {
      properties->data.UseZ = option_algorithms[properties->modeltype].Zdefault;
      sprintf(statusMessage,"Defaulting %s to %.*f\n",option_algorithms[properties->modeltype].UseZdes,properties->precision,option_algorithms[properties->modeltype].Zdefault);
      messageSet = 1;
    }

    if( option_algorithms[properties->modeltype].iUseB && properties->data.UseB == 0 &&
        option_algorithms[properties->modeltype].bBallow0Negative == 0 &&
        option_algorithms[properties->modeltype].Bdefault != 0 )
    {
      properties->data.UseB = option_algorithms[properties->modeltype].Bdefault;
      sprintf(statusMessage,"Defaulting %s to %.*f\n",option_algorithms[properties->modeltype].UseBdes,properties->precision,option_algorithms[properties->modeltype].Bdefault);
      messageSet = 1;
    }

    if( option_algorithms[properties->modeltype].iUseJ && properties->data.UseJ == 0 &&
        option_algorithms[properties->modeltype].bJallow0Negative == 0 &&
        option_algorithms[properties->modeltype].Jdefault != 0 )
    {
      properties->data.UseJ = option_algorithms[properties->modeltype].Jdefault;
      sprintf(statusMessage,"Defaulting %s to %.*f\n",option_algorithms[properties->modeltype].UseJdes,properties->precision,option_algorithms[properties->modeltype].Jdefault);
      messageSet = 1;
    }

    if(option_algorithms[properties->modeltype].iUseP && properties->data.UseP == 0 &&
       option_algorithms[properties->modeltype].bPallow0Negative == 0 &&
       option_algorithms[properties->modeltype].Pdefault != 0 )
    {
      properties->data.UseP = option_algorithms[properties->modeltype].Pdefault;
      sprintf(statusMessage,"Defaulting %s to %.*f\n",option_algorithms[properties->modeltype].UsePdes,properties->precision,option_algorithms[properties->modeltype].Pdefault);
      messageSet = 1;
    }

    if(option_algorithms[properties->modeltype].iUseQ && properties->data.UseQ == 0 &&
       option_algorithms[properties->modeltype].bQallow0Negative == 0 &&
       option_algorithms[properties->modeltype].Qdefault != 0 )
    {
      properties->data.UseQ = option_algorithms[properties->modeltype].Qdefault;
      sprintf(statusMessage,"Defaulting %s to %.*f\n",option_algorithms[properties->modeltype].UseQdes,properties->precision,option_algorithms[properties->modeltype].Qdefault);
      messageSet = 1;
    }

    if(option_algorithms[properties->modeltype].iUseR && properties->data.UseR == 0 &&
       option_algorithms[properties->modeltype].bRallow0Negative == 0 &&
       option_algorithms[properties->modeltype].Rdefault != 0 )
    {
      properties->data.UseR = option_algorithms[properties->modeltype].Rdefault;
      sprintf(statusMessage,"Defaulting %s to %.*f\n",option_algorithms[properties->modeltype].UseRdes,properties->precision,option_algorithms[properties->modeltype].Rdefault);
      messageSet = 1;
    }

    if(option_algorithms[properties->modeltype].iUseS && properties->data.UseS == 0 &&
       option_algorithms[properties->modeltype].bSallow0Negative == 0 &&
       option_algorithms[properties->modeltype].Sdefault != 0 )
    {
      properties->data.UseS = option_algorithms[properties->modeltype].Sdefault;
      sprintf(statusMessage,"Defaulting %s to %.*f\n",option_algorithms[properties->modeltype].UseSdes,properties->precision,option_algorithms[properties->modeltype].Sdefault);
      messageSet = 1;
    }

    if(option_algorithms[properties->modeltype].iUseT && properties->data.UseT == 0 &&
       option_algorithms[properties->modeltype].bTallow0Negative == 0 &&
       option_algorithms[properties->modeltype].Tdefault != 0 )
    {
      properties->data.UseT = option_algorithms[properties->modeltype].Tdefault;
      sprintf(statusMessage,"Defaulting %s to %.*f\n",option_algorithms[properties->modeltype].UseTdes,properties->precision,option_algorithms[properties->modeltype].Tdefault);
      messageSet = 1;
    }

    if( option_algorithms[properties->modeltype].bUsePound && properties->data.UsePound == 0 &&
        option_algorithms[properties->modeltype].bPoundallow0Negative == 0 &&
        option_algorithms[properties->modeltype].Pounddefault != 0)
    {
      properties->data.UsePound = option_algorithms[properties->modeltype].Pounddefault;

      if( option_algorithms[properties->modeltype].bUseStateNames == 0 )
      {
        sprintf(statusMessage,"Defaulting %s to %d\n",option_algorithms[properties->modeltype].UsePounddes,option_algorithms[properties->modeltype].Pounddefault);
        messageSet = 1;

      } else
      {
        sprintf(statusMessage,"Defaulting %s to %s\n",option_algorithms[properties->modeltype].UsePounddes,(char *)&(option_algorithms[properties->modeltype].StateNames[properties->data.UsePound-1]));
        messageSet = 1;
      }
    }

    /************************************************************************************/

    if( option_algorithms[properties->modeltype].iUseZ && option_algorithms[properties->modeltype].Zmax )
    {
      if( properties->data.UseZ > option_algorithms[properties->modeltype].Zmax )
      {
        properties->data.UseZ = option_algorithms[properties->modeltype].Zdefault;
      }
    }
    
    if( option_algorithms[properties->modeltype].iUseB && option_algorithms[properties->modeltype].Bmax )
    {
      if( properties->data.UseB > option_algorithms[properties->modeltype].Bmax )
      {
        properties->data.UseB = option_algorithms[properties->modeltype].Bdefault;
      }
    }

    if( option_algorithms[properties->modeltype].iUseJ && option_algorithms[properties->modeltype].Jmax )
    {
      if( properties->data.UseJ > option_algorithms[properties->modeltype].Jmax )
      {
        properties->data.UseJ = option_algorithms[properties->modeltype].Jdefault;
      }
    }

    if( option_algorithms[properties->modeltype].iUseP && option_algorithms[properties->modeltype].Pmax )
    {
      if( properties->data.UseP > option_algorithms[properties->modeltype].Pmax )
      {
        properties->data.UseP = option_algorithms[properties->modeltype].Pdefault;
      }
    }

    //

    if( option_algorithms[properties->modeltype].iUseQ && option_algorithms[properties->modeltype].Qmax )
    {
      if( properties->data.UseQ > option_algorithms[properties->modeltype].Qmax )
      {
        properties->data.UseQ = option_algorithms[properties->modeltype].Qdefault;
      }
    }

    if( option_algorithms[properties->modeltype].iUseR && option_algorithms[properties->modeltype].Rmax )
    {
      if( properties->data.UseR > option_algorithms[properties->modeltype].Rmax )
      {
        properties->data.UseR = option_algorithms[properties->modeltype].Rdefault;
      }
    }

    if( option_algorithms[properties->modeltype].iUseS && option_algorithms[properties->modeltype].Smax )
    {
      if( properties->data.UseS > option_algorithms[properties->modeltype].Smax )
      {
        properties->data.UseS = option_algorithms[properties->modeltype].Sdefault;
      }
    }

    if( option_algorithms[properties->modeltype].iUseT && option_algorithms[properties->modeltype].Tmax )
    {
      if( properties->data.UseT > option_algorithms[properties->modeltype].Tmax )
      {
        properties->data.UseT = option_algorithms[properties->modeltype].Tdefault;
      }
    }

    if( option_algorithms[properties->modeltype].bUsePound && option_algorithms[properties->modeltype].Poundmax )
    {
      if( properties->data.UsePound > option_algorithms[properties->modeltype].Poundmax )
      {
        properties->data.UsePound = option_algorithms[properties->modeltype].Pounddefault;

        if( option_algorithms[properties->modeltype].bUseStateNames == 1 ) 
        {
          //sprintf(statusMessage,"%s = %s",option_algorithms[properties->modeltype].UsePounddes,&(option_algorithms[properties->modeltype].StateNames[properties->data.UsePound-1]));

        } else
        {
          //sprintf(statusMessage,"%s = %d",option_algorithms[properties->modeltype].UsePounddes,properties->data.UsePound);
        }
        //messageSet = 1;
      }
    }

    /*******************************************************************/
    if(properties->format != DEMO_FUTURES && properties->format != DEMO_OPTIONS)
    {
      if( option_algorithms[properties->modeltype].assetClass == FUTURES_CLASS && 
          properties->format != CALENDAR_FUTURE &&
          properties->format != DECIMAL_FUTURE )
      {
        properties->format = CALENDAR_FUTURE;
        properties->decimalorcalendar=CALENDAR;
      }
    
      if( option_algorithms[properties->modeltype].assetClass != FUTURES_CLASS && 
         ( properties->format == CALENDAR_FUTURE ||
           properties->format == DECIMAL_FUTURE ))
      {
        properties->format = DECIMAL_OPTIONS;
        properties->decimalorcalendar = DECIMALS;
      }
    }

    if( properties->data.t[1] <= 0 )
    {
      properties->data.t[1] = 0;
      sprintf(statusMessage,"Please set second time variable to positive amount for pricing");
      messageSet = 1;
    }

    if( properties->data.t[2] <= 0 )
    {
      properties->data.t[2] = 0;
      sprintf(statusMessage,"Please set third time variable to positive amount for pricing");
      messageSet = 1;
    }

    if( properties->format != DEMO_OPTIONS )
    {
      if( option_algorithms[properties->modeltype].supportTime2 == 1 &&
          properties->format != DECIMAL_GREEKS &&
          properties->format != DECIMAL_OPTIONS )
      {
        properties->format = DECIMAL_OPTIONS;
        properties->decimalorcalendar = DECIMALS;
      }
    }

    if( option_algorithms[properties->modeltype].supportTime2 == 1 && 
        option_algorithms[properties->modeltype].perpetual == 0  &&
        properties->data.t[1] >= properties->data.t[0] )
    {
      properties->data.t[0] = properties->data.t[1]*2;
      sprintf(statusMessage,"Setting %s %.*f greater than %s %.*f",
              option_algorithms[properties->modeltype].supportTime1des,
              properties->precision,
              properties->data.t[0],
              option_algorithms[properties->modeltype].supportTime2des,
              properties->precision,
              properties->data.t[1]);
      messageSet = 1;
    }

    if( option_algorithms[properties->modeltype].supportTime3 == 1 && 
        properties->data.t[1] >= properties->data.t[2] )
    {
      properties->data.t[2] = properties->data.t[1]*2;
      sprintf(statusMessage,"Setting %s %.*f greater than %s %.*f",
              option_algorithms[properties->modeltype].supportTime3des,
              properties->precision,
              properties->data.t[2],
              option_algorithms[properties->modeltype].supportTime2des,
              properties->precision,
              properties->data.t[1]);
      messageSet = 1;
    }

    if( option_algorithms[properties->modeltype].supportVolatility == 1 && 
        properties->data.volatility < 0 )
    {
      properties->data.volatility = .25;
      sprintf(statusMessage,"Not set. Standard deviation forced to : %.3f",properties->data.volatility);
      messageSet = 1;
    }

    if( properties->decimalorcalendar == DECIMALS && properties->data.t[0] <= 0.0 )
    {
      //properties->data.t[0] = .25;
      properties->data.t[0] = 0;
      properties->data.te = 0;
      //sprintf(statusMessage,"Not set. Decimal time forced to: %.3f",properties->data.t[0]);
      sprintf(statusMessage,"Please set time to positive amount for pricing");
      messageSet = 1;
    }

    if( properties->realTimeBleeding == 1 )
    {
        time(&properties->updatetime);

        properties->data.te  = difftime(properties->updatetime,properties->starttime) / (double) 31536000;
        properties->data.te2 = difftime(properties->updatetime,properties->starttime2) / (double) 31536000;
        properties->data.te3 = difftime(properties->updatetime,properties->starttime3) / (double) 31536000;

        properties->data.te4 = difftime(properties->updatetime,properties->starttime4) / (double) 31536000;

    } else
    {
      properties->data.te  = 0;
      properties->data.te2 = 0;
      properties->data.te3 = 0;
      properties->data.te4 = 0;
    }

    if( properties->decimalorcalendar == CALENDAR )
    {
      properties->data.te = 0;
    }

    if( properties->customstrike <= 0 )
    {
      properties->customstrike = properties->data.price;
    }

    if( option_algorithms[properties->modeltype].bUsePound &&
        option_algorithms[properties->modeltype].bPoundallow0Negative != 1 && 
        properties->data.UsePound <= 0 )
    {
      properties->data.UsePound = 1;
      sprintf(statusMessage,"%s set to 1",option_algorithms[properties->modeltype].UsePounddes);
      messageSet = 1;
    }
    
    if( option_algorithms[properties->modeltype].supportDividend == 2 &&
        properties->format == DEMO_OPTIONS )
    {
      if(properties->data.dividend == 0)
      {
        properties->data.dividend = option_algorithms[properties->modeltype].defaultDividend;
      }
    }

    return messageSet;
}

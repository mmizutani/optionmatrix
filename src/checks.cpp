/* optionmatrix:                                                            */
/*                                                                          */
/* Options & Futures Matrix Modeler                                         */
/* View and Control Theoretical Option Chains                               */
/*                                                                          */
/* File: checks.cpp of optionmatrix                                         */
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
#include "main.h"
#include "license.h"

#ifdef HAVE_GETOPT_H
 #include <getopt.h>
#endif

int main(int argc, char **argv)
{
    setlocale(LC_ALL,"");
    programInits(&properties);
    props_defaults_options(&properties,1);

#ifdef HAVE_GETOPT_H

    process_arguments_checks(argc, argv, &properties.data.debug);

#endif

    program_source();
    program_check_pricing_models();

    exit(EXIT_SUCCESS);
}

#ifdef HAVE_GETOPT_H

void process_arguments_checks(int argc, char **argv, bool *debug)
{
    int c;
    int exit_program = FALSE;
    int exit_status = EXIT_FAILURE;

    *debug = false;

    static struct option long_options[] = {
        { "version", no_argument,       NULL, 'v' },
        { "help",    no_argument,       NULL, 'h' },
        { NULL,      no_argument,       NULL,  0  }
    };

    int option_index = 0;
    while ((c = getopt_long(argc, argv, "vh", long_options, &option_index)) != -1)
    {
        int this_option_optind = optind ? optind : 1;

        switch (c) {

        case 0:

            //printf ("option %s", long_options[option_index].name);

            if (optarg)
            {
                printf (" with arg %s not understood", optarg);

                exit_program = TRUE;
                exit_status = EXIT_FAILURE;
            }

            break;

         case 'v':

            program_version_checks();
            exit(EXIT_SUCCESS);

            break;

         case 'h':

            program_usage_checks();
            exit(EXIT_SUCCESS);

            break;

        case '?':

            exit_program = TRUE;
            exit_status = EXIT_FAILURE;

            break;

        default:

            printf("%s: ?? getopt returned character code 0%o ??\n", PACKAGE, c);

            exit_program = TRUE;
            exit_status = EXIT_FAILURE;

            break;

        }
    }
 
    if( optind < argc )
    {
        printf("\n%s: non-option argv-elements not understood: ", PACKAGE);

        while (optind < argc)
            printf ("%s ", argv[optind++]);
        printf ("\n");

        exit_program = TRUE;
        exit_status = EXIT_FAILURE;
    }

    if( exit_program )
     exit(exit_status);
    
    return;
}

void program_version_checks()
{

  printf("checks for gtk%s %s\n\n%s\n", PACKAGE, PACKAGE_VERSION, license2);

}

void program_usage_checks()
{

  printf("Usage: checks [OPTION] ...\n");

  printf("Confirm gtkoptionmatrix access to models source code.  Iterate through and compute all gtkoptionmatrix models with various strikes and inputs for testing.  Report total number of tests run.");

  printf("\n\n");
  printf("  -h, --help            display this help and exit\n");
  printf("  -v, --version         output version information\n");

  printf("\nReport %s bugs to %s\n", PACKAGE, PACKAGE_BUGREPORT);

}

#endif

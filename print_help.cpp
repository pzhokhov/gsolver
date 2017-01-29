#include "solver.h"

void print_helpmsg()
{
	printf(  "\nThis program solves generalized non-linear Shrodinger equation.");
	printf(  "\nCalling format:");
	printf("\n\nGNLSsolver scond nonlin [dump_prefix]");

	printf("\n\nscond       - file with starting condition (see help for startcondition)");

	printf("\n\nnonlin      - file with nonlinearity constants such as nonlinear refractive index");
	printf(  "\n              Raman constatns, ionization potential etc. See help for create_nl");
	printf(  "\ndump_prefix - prefix of the dump file. Dump names are formed as follows:");
	printf(  "\n              [dump_prefix]_id#######_r#");
	printf(  "\n              where id is followed by random integer number chosen at the beginning");
	printf(  "\n              r is followed by rank of the process ");
	printf(  "\n              Default dump prefix is GNLS_piece");

	printf(  "\n");
}


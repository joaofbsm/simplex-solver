/* Simplex Solver
 * Developed by Joao Francisco B. S. Martins <joaofbsm@dcc.ufmg.br>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "lalgebra.h"

/* PROGRAM FLOW *
 *
 * 1- Format LP to SEF
 * 2- Add operation register submatrix
 * 3- Create Auxiliar LP
 * 4- Solve Auxiliar LP
 * 5- If Auxiliar LP optimal value is:
 *    5.1- < 0, the LP is infeasible, so we should save the optimal solution for the dual of the 
 *         auxiliar LP(vector in [0->1][0->m - 2]) which is the infeasibility certificate
 *    5.2- = 0, the LP is feasible and so we need to extract the starting base from the resulting vector
 * 6- Now that we have a base to start with, we run the respective simplex. If mode is 1, we run the primal simplex. 
 *    If mode is 2, we run the chosen simplex. 
 *    6.1- Not regarding the chosen mode, we need to check in the simplex if we cannot proceed due to a column being 
 *         all negative in the primal simplex or a row being all positive in the dual simplex, the LP is unbounded. Using 
 *         previously created mechanisms you should be able to create an algorithm that finds the certificate.
 * 7- If we get to the end of the simplexes(canonical form, c and b >= 0) we can show the optimal solution and its value
 *    which is the first element of the last column.
 *    7.1- If mode equals 1 we need to print the optimal solution for the dual as well
 *    7.2- If mode equals 2 we need to print the tableaus step-by-step
 */

/* TODO: 
 * - Make matrix operations in another module
 * - Create struct of type LP with double** matrix, int m, int n, int type(-1 infeasible, 0 optimal, 1 unbounded)
 * - Refactor code optimizing it and renaming variables and functions
 * - Merge both simplex algorithm into only one
 */

/* PROBLEMS 
 * - Once again, when should we use an auxiliar LP in the second mode? When b is negative and c is positive in the tableau, we need the Auxiliar LP to get a viable base to start with.
 * - The LP can't be solved by the simplex, if we run the auxiliar first, sometimes. We need to correct that. There is a TODO explaining how in the format_canonical function
 * - opt_lista3ex1b.txt not working properly
 */

/* QUESTIONS
 * - Should our program adapt to solve, for instance, a problem with dual simplex where primal would be the right choice?
 */

int main(int argc, char* argv[]) {
    double** lp; // LP
    double** auxiliar_lp;
	int* base; // Bases are column numbers ordered by rows. If a column contains the base for the first restriction(first row of A), it is going to be on the first index of base and so forth
    int i, m, n, auxiliar_n, mode, simplex_result; // Matrix dimensions  
    char simplex_type; // Primal or Dual Simplex

    i = 0; // REMOVE THIS IN THE FINAL VERSION

    FILE* input; // Input file
    if(argc >= 2) { // Input file name has been given
		 input = fopen(argv[1], "r+");
    }
    else { // Default input file name
    	input = fopen("input.txt", "r+");
    }

    fscanf(input, "%d %d ", &m, &n);
    lp = allocate_matrix(m, n); // Allocate memory for matrix of dimensions m x n
    parse_input(input, lp, m, n);

    // Adds the slack variables for the problem, by formating it to the standard equalities form   
    lp = format_sef(lp, m, n);
    n += m - 1;

    format_tableau(lp, m, n);

//    printf("Standard Equality Form + Tableau(n = %d):\n", n);
//    print_matrix(lp, m, n);

   	// Now add the register of operations
	lp = add_operations_register(lp, m, n);
	n += m - 1; // Sums the new rows added in the previous operation to the total number of rows

//	printf("Operations Register(n = %d):\n", n);
//    print_matrix(lp, m, n);

	// User chooses the modus operandi
    printf("Escolha o modo de saída: ");
    scanf(" %d", &mode);
//    mode = 1;
//    printf("1\n\n");

    base = malloc((m - 1) * sizeof(int));

    switch(mode) {
	    case 1:

	    auxiliar_lp = create_auxiliar_lp(lp, m, n);
	    auxiliar_n = n + m - 1; // Value of n for the auxiliar PL with the operation register matrix on its side

//	    printf("\nAuxiliar(auxiliar_n = %d):\n", auxiliar_n);
//	    print_matrix(auxiliar_lp, m, auxiliar_n);

	    // Set the initial base for the auxiliar LP. j is the first row for the inserted columns
	    set_initial_base(auxiliar_lp, m, auxiliar_n, base);
//	    printf("Auxiliar Base:\n");
//	    for(i = 0; i < (m - 1); i++) {
//	    	printf("%d ", base[i]);
//	    }
//	    printf("\n\n");

	    // The base now is the final base of the auxiliar LP, which is a good one to begin the simplex with
		primal_simplex(auxiliar_lp, m, auxiliar_n, base, 0);
//		printf("Optimal Auxiliar LP:\n");
//		print_matrix(auxiliar_lp, m, auxiliar_n);

	    if(auxiliar_lp[0][auxiliar_n - 1] < 0) { // LP is infeasible
	    	printf("PL inviável, aqui está um certificado ");
	    	// The optimal solution for the dual of the auxiliar LP is a certificate of infeasibility for the original.
	    	print_output_vector(get_dual_optimal_solution(auxiliar_lp, m), m - 1);
	    	printf("\n");
	    }
	    else {
	    	simplex_result = primal_simplex(lp, m, n, base, 0);
	    	if(simplex_result > 0) { // LP is unbounded
	    		printf("PL ilimitada, aqui está um certificado ");
	    		print_output_vector(generate_unboundedness_certificate(lp, m, n, simplex_result, base), (n - 1 - (m - 1) - (m - 1)));
	    		printf("\n");
	    	} 
	    	else { // LP is optimal
	    		printf("Solução ótima x = ");
	    		print_output_vector(get_primal_optimal_solution(lp, m, n, base), (n - 1 - (m - 1) - (m - 1)));
	    		printf(", com valor objetivo %g, e solução dual y = ", round(lp[0][n - 1] * 100000) / 100000);
	    		print_output_vector(get_dual_optimal_solution(lp, m), m - 1);
	    		printf("\n");
	    	}
	    }

		free(auxiliar_lp);

	    break;

	    case 2:
	    	printf("Você gostaria de resolver pelo simplex (P)rimal ou (D)ual? ");
	    	scanf(" %c", &simplex_type);

	    	// Set bases to the slack variables
			set_initial_base(lp, m, n, base);

	    	switch(simplex_type) {
	    		case 'P':
	    		simplex_result = primal_simplex(lp, m, n, base, 1);
	    		break;

	    		case 'D':
	    		simplex_result = dual_simplex(lp, m, n, base, 1);
	    		break;

	    		default:
				printf("Erro: Opção Inválida.\n");
	    	}
	    break;

		default:
		printf("Erro: Opção Inválida.\n");
	}

	fclose(input);

	free(lp);
	free(base);

    return 0;
}

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
 * - Create struct of type LP with double** matrix, int m, int n, int type(-1 infeasible, 0 optimal, 1 unbounded)
 */

/* PROBLEMS 
 * - Can't adapt Primal -> Dual
 */

/* QUESTIONS
 * - Program adapts itself to solve LP's if the wrong simplex is requested. It works for Primal Simplex but 
 *   not for Dual(The chosen bases on the Auxiliar LP are not good enough). In those cases, the solution 
 *   starts at the optimal but walks back to the original form -> ASK TEACHER FOR HELP
 * - Also check if the program should do this adaptation. For instance, the answer for opt_lista3ex1b.txt
 *   is tableau sequence is different for the dual input with this implementation
 * - Should find_non_zero_element be > 0?
 * - What else needs to be printed with the code?
 */

int main(int argc, char* argv[]) {
    // Linear Programming represented as a matrix similar to the tableau
    double** lp; 

    // Auxiliar LP built upon the original LP
    double** auxiliar_lp;

    // Bases are column numbers ordered by rows. If a column contains the base for the first 
    // restriction(first row of A), it is going to be on the first index of base and so forth
	int* base; 

	// m and n are the dimensions of the LP. auxiliar_n is the columns dimension of the auxiliar_lp. 
	// mode is the mode chosen by the user. simplex_result is the return value of the simplex algorithms
    int m, n, auxiliar_n, mode, simplex_result;

    // User choice for primal or dual simplex in mode 2
    char simplex_type; 

    // Input file
    FILE* input; 
    if(argc >= 2) { // Input file name has been given
		 input = fopen(argv[1], "r+");
    }
    else { // Default input file name
    	input = fopen("input.txt", "r+");
    }

    fscanf(input, "%d %d ", &m, &n); // Reads LP dimensions
    lp = allocate_matrix(m, n); // Allocate memory for matrix of dimensions m x n
    parse_input(input, lp, m, n); // Fill the allocated matrix with the input

    lp = format_sef(lp, m, n); // Adds the slack variables for the problem by formating it to the standard equalities form   
    n += m - 1; // (m - 1) columns were added to the matrix

    format_tableau(lp, m, n); // Negates the first row for the tableau

	lp = add_operations_register(lp, m, n); // Adds the operation register matrix to the left of the LP
	n += m - 1;

	// User chooses the modus operandi
    printf("Escolha o modo de saída: ");
    scanf(" %d", &mode);

    auxiliar_lp = NULL;
    base = malloc((m - 1) * sizeof(int)); // Base will always be a vector with (m - 1) columns because this is the rank of the matrix

    switch(mode) {
	    case 1:
		    auxiliar_lp = create_auxiliar_lp(lp, m, n); 
		    auxiliar_n = n + m - 1; // Auxiliar LP creates (m - 1) new columns in A

		    set_initial_base(auxiliar_lp, m, auxiliar_n, base); // Set the initial base for the auxiliar LP

			primal_simplex(auxiliar_lp, m, auxiliar_n, base, 0); // Runs simplex for Auxiliar LP but doesn't print the output

		    if(auxiliar_lp[0][auxiliar_n - 1] < 0) { // LP is infeasible
		    	printf("PL inviável, aqui está um certificado ");
		    	// The optimal solution for the dual of the auxiliar LP is a certificate of infeasibility for the original LP
		    	print_output_vector(get_dual_optimal_solution(auxiliar_lp, m), m - 1);
		    	printf("\n");
		    }
		    else {
		    	// The base now is the final base of the auxiliar LP, which is a good one to begin the simplex with
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
	    break;

	    case 2:
	    	printf("Você gostaria de resolver pelo simplex (P)rimal ou (D)ual? ");
	    	scanf(" %c", &simplex_type);

	    	switch(simplex_type) {
	    		case 'P':
	    			// If b has some negative entry, use auxiliar LP to find a good base of columns to start the simplex with
		    		if(is_b_negative(lp, m, n)) { 
		    			auxiliar_lp = create_auxiliar_lp(lp, m, n);
		    			auxiliar_n = n + m - 1;
		    			set_initial_base(auxiliar_lp, m, auxiliar_n, base);
		    			primal_simplex(auxiliar_lp, m, auxiliar_n, base, 0);
		    		}
		    		else {
		    			// Set base columns to the slack variables
						set_initial_base(lp, m, n, base);
		    		}

		    		simplex_result = primal_simplex(lp, m, n, base, 1);
	    		break;

	    		case 'D':
	    			// The same procedure we did for primal simplex, but this time it is done with the vector c instead of b
		    		if(is_c_negative(lp, m, n)) {
		    			auxiliar_lp = create_auxiliar_lp(lp, m, n);
		    			auxiliar_n = n + m - 1;
		    			set_initial_base(auxiliar_lp, m, auxiliar_n, base);
		    			primal_simplex(auxiliar_lp, m, auxiliar_n, base, 0);
		    		}
		    		else {
		    			// Set base columns to the slack variables
						set_initial_base(lp, m, n, base);
		    		}

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
	if(auxiliar_lp != NULL) { // If it was used to solve the LP, we need to free it
		free(auxiliar_lp); 
	}

    return 0;
}
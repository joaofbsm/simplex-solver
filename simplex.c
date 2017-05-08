/* Simplex Solver
 * Developed by Joao Francisco B. S. Martins <joaofbsm@dcc.ufmg.br>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "lalgebra.h"

int main(int argc, char* argv[]) {
  // Linear Programming represented as a matrix similar to the tableau
  double** lp; 

  // Auxiliar LP built upon the original LP
  double** auxiliar_lp;

  // Bases are column numbers ordered by rows. If a column contains the base for the first 
  // restriction(first row of A), it is going to be on the first index of base and so forth
  int* base; 

  // m and n are the dimensions of the LP. auxiliar_n is the columns dimension for the auxiliar_lp. 
  // mode is the mode chosen by the user. simplex_result is the return value of the simplex algorithms
  int m, n, auxiliar_n, mode, simplex_result;

  // User choice for primal or dual simplex in mode 2
  char simplex_type; 

  // Receive and discard the string "modo". We make this to facilitate the parsing
  char string_modo[4];

  // Input file
  FILE* input; 
  if(argc >= 2) { // Input file name has been given
     input = fopen(argv[1], "r+");
  }
  else { // Default input file name
    input = fopen("input.txt", "r+");
  }

  // Gets the modus operandi
  fscanf(input, "%s %d", string_modo, &mode);

  // Primal or dual simplex
  if(mode == 2) { 
    fscanf(input, " %c", &simplex_type); 
  }

  fscanf(input, "%d %d ", &m, &n); // Reads LP dimensions
  m += 1; // Adjust m so it corresponds to the tableau dimensions
  n += 1; // Adjust n so it corresponds to the tableau dimensions
  lp = allocate_matrix(m, n); // Allocate memory for matrix of dimensions m x n
  parse_input(input, lp, m, n); // Fill the allocated matrix with the input

  lp = format_sef(lp, m, n); // Adds the slack variables for the problem by formating it to the standard equalities form   
  n += m - 1; // (m - 1) columns were added to the matrix

  format_tableau(lp, m, n); // Negates the first row for the tableau

  lp = add_operations_register(lp, m, n); // Adds the operation register matrix to the left of the LP
  n += m - 1;

  auxiliar_lp = NULL;
  // Base will always be a vector with (m - 1) columns because this is the rank of the matrix
  base = malloc((m - 1) * sizeof(int)); 

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
          // We can only run the dual simplex if we have a positive c vector in the tableau.
          // That is not entirely true, but for the scope of this assignment it will
          if(is_c_positive(lp, m, n)) {
            set_initial_base(lp, m, n, base);
            simplex_result = dual_simplex(lp, m, n, base, 1);
          }
          else {
            printf("Não foi possível rodar o simplex dual com a PL dada.\n");
          }
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
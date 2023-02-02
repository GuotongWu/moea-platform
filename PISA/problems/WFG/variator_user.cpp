/*========================================================================
  PISA  (www.tik.ee.ethz.ch/pisa/)
 
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich
 
  ========================================================================
  PISALIB 
  
  Functions the user must implement
  
  C file.
  
  file: variator_user.c
  author: Fabian Landis, flandis@ee.ethz.ch

  revision by: Stefan Bleuler, bleuler@tik.ee.ethz.ch
  last change: $date$
  
  ========================================================================
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "wfg_interface.h"
#include "variator.h"
#include "variator_user.h"



/*----------------------| global variables |---------------------------*/

/* declared extern in variator_user.h and used in other files as well */

char paramfile[FILE_NAME_LENGTH]; /* file with local parameters */

char *log_file = "variator_error.log";
/**** Change the definition of 'log_file' if you want another name! Or
      set to NULL if you want to write error messages to stderr. */

/* local parameters from paramfile*/
char problem[FILE_NAME_LENGTH]; 
int seed;   /* seed for random number generator */
int number_decision_variables; /* length of the binary string */
int number_position_variables;
int number_distance_variables;
int maxgen; /* maximum number of generations (stop criterion) */
int gen;
char outfile[FILE_NAME_LENGTH]; /* output file for last population */
double individual_mutation_probability;
double individual_recombination_probability;
double variable_mutation_probability;
double variable_swap_probability;
double variable_recombination_probability;
double eta_mutation;
double eta_recombination;
int use_symmetric_recombination;


/*-------------------------| individual |-------------------------------*/

void free_individual(individual *ind) 
/* Frees the memory for one indiviual.

   post: memory for ind is freed
*/
{
	/**********| added for WFG |**************/
	if (ind == NULL)
		return;

	free(ind->x);
	free(ind->f);

	/**********| addition for WFG end |*******/

	free(ind);
}

double get_objective_value(int identity, int i)
/* Gets objective value of an individual.

   pre: 0 <= i <= dimension - 1 (dimension is the number of objectives)

   post: Return value == the objective value number 'i' in individual '*ind'.
         If no individual with ID 'identity' return value == -1. 
*/   
{
	double objective_value = -1.0;

	assert(0 <= i && i < dimension); /* asserting the pre-condition */

	/**********| added for DTLZ |**************/
	individual *temp_ind;

	if (i < 0 || i > (dimension - 1))
		return(-1);

	temp_ind = get_individual(identity);
	if (temp_ind == NULL)
		return(-1);

	objective_value = temp_ind->f[i];
	/**********| addition for WFG end |*******/

	return (objective_value);
}

/*-------------------------| statemachine functions |-------------------*/

int state0() 
/* Do what needs to be done in state 0.

   pre: The global variable 'paramfile' contains the name of the
        parameter file specified on the commandline.
        The global variable 'alpha' contains the number of indiviuals
        you need to generate for the initial population.
                
   post: Optionally read parameter specific for the module.
         Optionally do some initialization.
         Initial population created.
         Information about initial population written to the ini file
         using write_ini().
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
    
     int result; /* stores return values of called functions */
     int *initial_population; /* storing the IDs of the individuals */
     
     /**********| added for WFG |**************/
     int i;
     /**********| addition for WFG end |*******/
     
     initial_population = (int *) malloc(alpha * sizeof(int)); 
     if (initial_population == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }
     
     /**********| added for WFG |**************/
     result = read_local_parameters();

     if (result != 0)
     { 
    	 log_to_file(log_file, __FILE__, __LINE__,
    			 "couldn't read local parameters");
    	 return (1);
     }

     /* initializing the first alpha individuals */
     for(i = 0; i < alpha; i++)
     {
    	 individual *ind = new_individual();
    	 eval(ind);
    	 initial_population[i] = add_individual(ind);
    	 if(initial_population[i] == -1)
    		 return(1);
     } 

     gen = 1;

     /**********| addition for WFG end |*******/

     result = write_ini(initial_population);
     if (result != 0)
     { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't write ini");
          free(initial_population);
          return (1);
     }

     free(initial_population);
     return (0);
}



int state2()
/* Do what needs to be done in state 2.

   pre: The global variable 'mu' contains the number of indiviuals
        you need to read using 'read_sel()'.
        The global variable 'lambda' contains the number of individuals
        you need to create by variation of the individuals specified the
        'sel' file.
        
   post: Optionally call read_arc() in order to delete old uncessary
         individuals from the global population.
         read_sel() called
         'lambda' children generated from the 'mu' parents
         Children added to the global population using add_individual().
         Information about children written to the 'var' file using
         write_var().
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     int *parent_identities, *offspring_identities; /* array of identities */
     int result; /* stores return values of called functions */

     parent_identities = (int *) malloc(mu * sizeof(int)); 
     if (parent_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }

     offspring_identities = (int *) malloc(lambda * sizeof(int)); 
     if (offspring_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }
     
     result = read_sel(parent_identities);
     if (result != 0) /* if some file reading error occurs, return 2 */
          return (2);

     result = read_arc(); 
     if (result != 0) /* if some file reading error occurs, return 2 */
          return (2);

     /**********| added for WFG |**************/

     result = variate(parent_identities, offspring_identities);
     if (result != 0)
    	 return (1);

     gen++;

     /**********| addition for WFG end |*******/

     result = write_var(offspring_identities);
     if (result != 0)
     { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't write var");
          free(offspring_identities);
          free(parent_identities);
          return (1);
     }

     free(offspring_identities);
     free(parent_identities);
     return (0);
}
 

int state4() 
/* Do what needs to be done in state 4.

   pre: State 4 means the variator has to terminate.

   post: Free all memory.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
	/**********| added for WFG |**************/

	int result;
	result = read_arc();

	if (0 == result) /* arc file correctly read
	                         this means it was not read before,
	                         e.g., in a reset. */
	{
		write_output_file();
	}

	/**********| addition for WFG end |*******/
     
     return (0);
}


int state7()
/* Do what needs to be done in state 7.

   pre: State 7 means that the selector has just terminated.

   post: You probably don't need to do anything, just return 0.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     return(0);  
}


int state8()
/* Do what needs to be done in state 8.

   pre: State 8 means that the variator needs to reset and get ready to
        start again in state 0.

   post: Get ready to start again in state 0. 
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
	/**********| added for DTLZ |**************/

	int result;
	gen = 1;
	result = read_arc();

	if (0 == result) /* arc file correctly read
	                       this means it was not read before */
	{
		write_output_file();
	}

	/**********| addition for DTLZ end |*******/

	return (0);
}


int state11()
/* Do what needs to be done in state 11.

   pre: State 11 means that the selector has just reset and is ready
        to start again in state 1.

   post: You probably don't need to do anything, just return 0.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     return (0);  
}


int is_finished()
/* Tests if ending criterion of your algorithm applies.

   post: return value == 1 if optimization should stop
         return value == 0 if optimization should continue
*/
{    
     /**** Here comes your code for checking the termination criteria. */
	/**********| added for WFG |**************/	
	return (gen >= maxgen);
	/**********| addition for WFG end |*******/
}

/**********| added for DTLZ |**************/

int read_local_parameters()
{
     FILE *fp;
     char str[CFG_NAME_LENGTH];

     /* reading parameter file with parameters for selection */
     fp = fopen(paramfile, "r"); 
     assert(fp != NULL);

     if(dimension < 0)
     {
          log_to_file(log_file, __FILE__, 
                      __LINE__, "can't handle that dimension");
          return(1);
     } 

     if(mu != lambda)
     {
          log_to_file(log_file, __FILE__, 
                      __LINE__, "can't handle mu != lambda");
          return(1);
     }

     fscanf(fp, "%s", str);
     assert(strcmp(str, "problem") == 0);
     fscanf(fp, "%s", problem); /* fscanf() returns EOF if
                                   reading failed. */

     fscanf(fp, "%s", str);
     assert(strcmp(str, "seed") == 0);
     fscanf(fp, "%d", &seed);
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "number_decision_variables") == 0);
     fscanf(fp, "%d", &number_decision_variables);
          
     fscanf(fp, "%s", str);
     assert(strcmp(str, "number_position_variables") == 0);
     fscanf(fp, "%d", &number_position_variables);
          
     fscanf(fp, "%s", str);
     assert(strcmp(str, "number_distance_variables") == 0);
     fscanf(fp, "%d", &number_distance_variables);
     
     if( number_decision_variables < dimension )
     	fprintf(stderr, "Number of decision variables %d smaller than number of dimensions", number_decision_variables, dimension );
     if( number_decision_variables != number_position_variables +
     	number_distance_variables )
     	fprintf(stderr, "Number of position variables (%d) + number of distance variables (%d) does not equal number_decision_variables (%d)\n",
     	number_position_variables, number_distance_variables, number_decision_variables );
     
     assert( number_decision_variables >= dimension );
     assert( number_position_variables >= 0 );
     assert( number_distance_variables >= 0 );
     assert( number_decision_variables == number_position_variables + 
    	 	number_distance_variables );
     
     

     fscanf(fp, "%s", str);
     assert(strcmp(str, "maxgen") == 0);
     fscanf(fp, "%d", &maxgen);
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "outputfile") == 0);
     fscanf(fp, "%s", outfile); /* fscanf() returns EOF if
                                   reading failed. */

     fscanf(fp, "%s", str);
     assert(strcmp(str, "individual_mutation_probability") == 0);
     fscanf(fp, "%le", &individual_mutation_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "individual_recombination_probability") == 0);
     fscanf(fp, "%le", &individual_recombination_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "variable_mutation_probability") == 0);
     fscanf(fp, "%le", &variable_mutation_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "variable_swap_probability") == 0);
     fscanf(fp, "%le", &variable_swap_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "variable_recombination_probability") == 0);
     fscanf(fp, "%le", &variable_recombination_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "eta_mutation") == 0);
     fscanf(fp, "%le", &eta_mutation);
     assert(eta_mutation >= 0);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "eta_recombination") == 0);
     fscanf(fp, "%le", &eta_recombination);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "use_symmetric_recombination") == 0);
     fscanf(fp, "%d", &use_symmetric_recombination);

     
     srand(seed); /* seeding random number generator */

     fclose(fp);

     return (0);
}

/* Generate a random integer. */
int irand(int range)
{
	int j;
	j=(int) ((double) range * (double) rand() / (RAND_MAX + 1.0));
	return (j);
}

/* Generate a random double. */
double drand(double range)
{
	double j;
	j=(range * (double) rand() / (RAND_MAX + 1.0));
	return (j);
}

/* Determines the objective value based on DTLZ */
int eval(individual *ind)
{
	double fit[ dimension ];
	wfg_eval( ind->x, number_decision_variables, number_position_variables,
			  dimension, problem, fit );
	for( int i = 0; i < dimension; i++ )
		ind->f[i] = fit[i]; 
	return 0;
}

/* Writes the index, objective values and bit string of
   all individuals in global_population to 'out_filename'. */
void write_output_file()
{
     int j, current_id;
     FILE *fp_out;
     individual *temp;
     
     fp_out = fopen(outfile, "w");
     assert(fp_out != NULL);

     current_id = get_first();

     while (current_id != -1)
     {       
	  temp = get_individual(current_id);
          fprintf(fp_out, "%d ", current_id); /* write index */
	  for (j = 0; j < dimension; j++)
	  {
	      fprintf(fp_out, "%f ", temp->f[j]);
	  }
          
	  for (j = 0; j < temp->n; j++)
          {
               fprintf(fp_out, "%f ", temp->x[j]);
          }
          fprintf(fp_out, "\n");
	  current_id = get_next(current_id);
     }

     fclose(fp_out);
}

/* create a random new individual and allocate memory for it,
   returns a pointer to the new individual */
individual *new_individual()
{
     individual *return_ind;
     int i;

     return_ind = (individual *) malloc(sizeof(individual));
     return_ind->x = (double *) malloc(sizeof(double) * number_decision_variables);
     return_ind->f = (double *) malloc(sizeof(double) * dimension);
 
     for (i = 0; i < number_decision_variables; i++)
     {
	 return_ind->x[i] = drand(1.0);
     }
     
     return_ind->n = number_decision_variables;

     for (i = 0; i < dimension; i++)
     {
	 return_ind->f[i] = 0;
     }
     
     return (return_ind);
}

/* Performs variation. */
int variate(int *selected, int *result_ids)
{
	int result, i, k;

	result = 1;

	/* copying all individuals from selected */
	for(i = 0; i < mu; i++)
	{
		result_ids[i] = 
			add_individual(copy_individual(get_individual(selected[i])));
		if(result_ids[i] == -1)
		{
			log_to_file(log_file, __FILE__, __LINE__,
					"copying + adding failed");
			return (1);
		}
	}

	/* if odd number of individuals, last one is
        left as is */
	if((((double)mu/2) - (int)(mu/2)) != 0) k = mu - 1; 
	else k = mu;

	/* do recombination */
	for(i = 0; i < k; i+= 2)
	{  
		if (drand(1) <= individual_recombination_probability)
		{
			if (variable_swap_probability > 0)
			{
				result = uniform_crossover
				(get_individual(result_ids[i]),
						get_individual(result_ids[i + 1]));
				if (result != 0)
					log_to_file(log_file, __FILE__, 
							__LINE__, "recombination failed!");
			}

			if (variable_recombination_probability > 0)
			{
				result = sbx
				(get_individual(result_ids[i]), 
						get_individual(result_ids[i + 1]));
				if (result != 0)
					log_to_file(log_file, __FILE__, 
							__LINE__, "recombination failed!");
			}
		}
	}

	/* do mutation */
	for(i = 0; i < mu; i++)
	{
		if (drand(1) <= individual_mutation_probability)
		{ 
			if (variable_mutation_probability > 0)
			{
				result = mutation(get_individual(result_ids[i]));
				if(result != 0)
					log_to_file(log_file, __FILE__, __LINE__,
							"mutation failed!");
			}
		}
	}

	/* do evaluation */
	for(i = 0; i < mu; i++)
	{
		int result;
		result = eval(get_individual(result_ids[i]));
	}

	return (0);
}

int mutation(individual *ind)
{
	int i;

	if (ind == NULL)
	{
		return (1);
	}

	for (i = 0; i < ind->n; i++)
	{
		if (drand(1) <= variable_mutation_probability)
		{
			double eta = eta_mutation;
			double u = drand(1.0);
			double delta = 0;
			double x = ind->x[i];
			double lb = 0;    /* lower bound of variable i */
			double ub = 1;    /* upper bound of variable i */
			double diff = ub - lb;  /* range of variable i */
			double maxmut0 = x - lb;
			double maxmut = ub - x;
			double delta_max = maxmut0 / diff;
			if (maxmut0 > maxmut)
			{
				delta_max = maxmut / diff;
			}

			if (u < 0.5)
			{
				double b =  2*u + (1-2*u)*(pow(1-delta_max,(eta+1)));
				delta = pow(b,(1.0/(eta+1))) - 1.0;
			}
			else
			{
				double b = 2*(1-u) + 2*(u-0.5)*(pow(1-delta_max,(eta+1)));
				delta = 1.0 - pow(b,(1.0/(eta+1)));
			}
			if (delta > delta_max)  /* machine accuracy problem */
				delta = delta_max;
			else if (delta < -delta_max)
				delta = -delta_max;

			ind->x[i] = x + delta * diff;
		}
	}

	return (0);
}

int uniform_crossover(individual *ind1, individual *ind2)
{
	int i;

	for (i = 0; i < ind2->n; i++)
	{
		if (drand(1) <= variable_swap_probability) /* switch variable */
		{
			double x = ind2->x[i];
			ind2->x[i] = ind1->x[i];
			ind1->x[i] = x;
		} 
	}  

	return (0);
}

int sbx(individual *ind1, individual *ind2)
{
	int i;

	for (i = 0; i < ind2->n; i++)
	{
		if (drand(1) <= variable_recombination_probability)  
		{
			double di = eta_recombination; /* distribution index */
			int bounded = 1;
			double lb = 0;    /* lower bound of variable i */
			double ub = 1;    /* upper bound of variable i */	     
			double u = drand(1);
			double b0=0, b1=0;   /* spread factors */
			double x0 = ind1->x[i];
			double x1 = ind2->x[i];
			double bl=0, bu=0, p_bl=0, p_bu=0, bll=0, buu=0, blll=0, buuu=0;
			double dx = 0;
			double u0=0, u1=0;

			/* calculate spread factor(s) b0, b1 */ 
			if (bounded == 1)
			{
				dx = fabs(x1-x0);   /* difference of x values */
				if (dx > 0)
				{
					bl = (x0 + x1 - 2*lb) / dx;
					bu = (2*ub - x0 - x1) / dx;
					bll = (x0 + x1 - 2*(x0-lb)) / dx;
					buu = (2*(ub-x1)-x0-x1) / dx;
					if (x0 < x1)
					{
						blll = 1 + 2 * (x0 - lb) / dx;
						buuu = 1 + 2 * (ub - x1) / dx;
					}
					else
					{
						blll = 1 + 2 * (x1 - lb) / dx;
						buuu = 1 + 2 * (ub-x0) / dx;
					}

					bl = blll; /* take Deb's version (numerically stable) */
					bu = buuu;

					/* switch off symmetric recombination to avoid
					 * getting stuck on a line where one parameter
					 * equals an extreme value. */
					 if (use_symmetric_recombination)
					 {
						 if (bl < bu)  /* symmetric distribution, like Deb */
							 bu = bl;
						 else
							 bl = bu;
						 assert (b0 == b1);
					 }
					 assert(bl > 0 && bu > 0);
					 p_bl = 1 - 1/(2*pow(bl,di+1));
					 p_bu = 1 - 1/(2*pow(bu,di+1));
				}
				else
				{
					p_bl = 1;
					p_bu = 1;
				}
				u0 = u*p_bl;
				u1 = u*p_bu;
				if (u0<=0.5)
					b0 = pow(2*u0,1/(di+1));
				else
					b0 = pow(0.5/(1-u0),1/(di+1));
				if (u1<=0.5)
					b1 = pow(2*u1,1/(di+1));
				else
					b1 = pow(0.5/(1-u1),1/(di+1));
				assert(dx==0 || (b0<=bl && b1<=bu)); /* machine accuracy */
			}
			else
			{
				if (u<=0.5)
					b0 = pow(2*u,1/(di+1));
				else
					b0 = pow(0.5/(1-u),1/(di+1));
				b1 = b0;
			}

			if (x0 < x1)
			{
				ind1->x[i] = 0.5*(x0+x1 + b0*(x0-x1));
				ind2->x[i] = 0.5*(x0+x1 + b1*(x1-x0));
			}
			else
			{
				ind1->x[i] = 0.5*(x0+x1 + b1*(x0-x1));
				ind2->x[i] = 0.5*(x0+x1 + b0*(x1-x0));
			}
		}
	}  

	return (0);
}

/* copy an individual and return the pointer to it */
individual *copy_individual(individual *ind)
{
     individual *return_ind;
     int i;

     return_ind = (individual *) malloc(sizeof(individual));
     return_ind->x = (double *)malloc(sizeof(double) * number_decision_variables);
     return_ind->f = (double *) malloc(sizeof(double) * dimension);
     
     for (i = 0; i < number_decision_variables; i++)
          return_ind->x[i] = ind->x[i];

     for (i = 0; i < dimension; i++)
	  return_ind->f[i] = ind->f[i];

     return_ind->n = ind->n;

     return(return_ind);
}

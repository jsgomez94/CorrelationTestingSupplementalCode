#################################################
#################################################
#################################################
##
## Important function: CreatingParameters()
##
#################################################
#################################################
#################################################


#################################################
#################################################
## CreateParameters:
##    Function that, given a particular simulation
##    ID id_task, returns the simulation parameters
##    corresponding to such task ID. 
##
##  INPUTS
##    id_task  : ID of the current simulation task.
##                If id_task==0, is a small debugging example.
##                If id_task%in$1:288, simulations with
##                  different choices of p,T,n, ph, pnh, etc.
##    
##  OUTPUT:
##
create_parameters <- function(id_task) {
  ## id_task = 0 
  ##  corresponds to a reduced experiment that is useful
  ##  for debugging the code.
  if(id_task == 0) {
    args_list <- list(
      nnn             = 50,
      model           = 2,
      sp              = 0,
      p               = 10,
      id_task         = 0)

    return(args_list)
  }
  ## id_task in 1-288
  ##  corresponds to the parameters used for our 
  ##  systematic simulations.

  ## TABLE OF ALL PARAMETER COMBINATIONS.
  sim_par_table <- expand.grid(
    nnn             = c(40, 80, 160, 320),
    model           = c(1, 2, 3, 4, 5),
    sp              = c(0, 1, 2),
    p               = c(15, 30, 60, 120, 240, 480, 960))

  ## Function returns row of index id_task.
  args <- sim_par_table[id_task, ]
  args_list <- list()
  for (i in 1:ncol(args)) {
    args_list[[i]] <- args[1,i]
  }
  names(args_list) <- colnames(args)
  
  args_list$id_task <- id_task
  
  return(args_list)
}


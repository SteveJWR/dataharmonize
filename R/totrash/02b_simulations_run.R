library(parallel)
numCores <- detectCores()
numCores
source("R/01_functions.R")
source("R/02a_simulations_setup.R")


# setting the seed across simulations
RNGkind("L'Ecuyer-CMRG")
set.seed(1)






RUN_PARALLEL = FALSE
timeout = 5000





##### ----- intrinsic variability simulations model 1 -----

intrinsic.variability.model1.simulation.results.array <- array(NA, dim = c(n.sims, length(h.set), length(ker.set), length(mu.set)))


intrinsic.variability.model1.simulation.results <- simulate_experiment(RUN_PARALLEL,
                                                                       simulation_intrinsic_variability_model1,
                                                                       n.sims)


for(i in 1:n.sims){
  for(q in 1:nrow(hyper.param.idx)){
    j = hyper.param.idx[q,1]
    k = hyper.param.idx[q,2]
    l = hyper.param.idx[q,3]
    intrinsic.variability.model1.simulation.results.array[i,j,k,l] <- intrinsic.variability.model1.simulation.results[[i]][q]
  }
}

saveRDS(intrinsic.variability.model1.simulation.results.array,paste0("Data/results/intrinsic_variability_simulation_model1.rds"))
print("intrinsic variability simulations model 1 complete")




##### ----- intrinsic variability simulations model 2 -----

intrinsic.variability.model2.simulation.results.array <- array(NA, dim = c(n.sims, length(h.set), length(ker.set), length(mu.set)))


intrinsic.variability.model2.simulation.results <- simulate_experiment(RUN_PARALLEL,
                                                                       simulation_intrinsic_variability_model2,
                                                                       n.sims)


for(i in 1:n.sims){
  for(q in 1:nrow(hyper.param.idx)){
    j = hyper.param.idx[q,1]
    k = hyper.param.idx[q,2]
    l = hyper.param.idx[q,3]
    intrinsic.variability.model2.simulation.results.array[i,j,k,l] <- intrinsic.variability.model2.simulation.results[[i]][q]
  }
}



saveRDS(intrinsic.variability.model2.simulation.results.array,paste0("Data/results/intrinsic_variability_simulation_model2.rds"))

print("intrinsic variability simulations model 2 complete")





## -------- mu selection --------
## --------    Model 1   --------



smoothing.selection.model1.simulation.results.array <- array(NA, dim = c(n.sims, length(mu.set.long)))

smoothing.selection.model1.simulation.results <- simulate_experiment(RUN_PARALLEL,
                                                                     simulation_smoothing_selection_model1,
                                                                     n.sims)

for(i in 1:n.sims){
  for(q in 1:length(mu.set.long)){
    smoothing.selection.model1.simulation.results.array[i,q] <- smoothing.selection.model1.simulation.results[[i]][q]
  }
}



saveRDS(smoothing.selection.model1.simulation.results.array,paste0("Data/results/smoothing_selection_simulation_model1.rds"))


print("mu selection simulations model 1 complete")


## --------    Model 2   --------


smoothing.selection.model2.simulation.results.array <- array(NA, dim = c(n.sims, length(mu.set.long)))

smoothing.selection.model2.simulation.results <- simulate_experiment(RUN_PARALLEL,
                                                                     simulation_smoothing_selection_model2,
                                                                     n.sims)

for(i in 1:n.sims){
  for(q in 1:length(mu.set.long)){
    smoothing.selection.model2.simulation.results.array[i,q] <- smoothing.selection.model2.simulation.results[[i]][q]
  }
}

saveRDS(smoothing.selection.model2.simulation.results.array,paste0("Data/results/smoothing_selection_simulation_model2.rds"))


print("mu selection simulations model 2 complete")





## -------- Conversion Cross Entropy --------



conversion.ce.simulation.results.array <- array(NA, dim = c(n.sims.conversion, length(h.set.conversion), length(ker.set), length(mu.set.conversion)))
conversion.ce.simulation.ml.results.array <- array(NA, dim = c(n.sims.conversion, length(h.set.conversion), length(ker.set)))


conversion.ce.simulation.results <- simulate_experiment(RUN_PARALLEL,
                                                        simulation_conversion_cross_entropy,
                                                        n.sims.conversion)
saveRDS(conversion.ce.simulation.results, "Data/results/conversion_ce_simulations_list1.rds")



for(i in 1:n.sims.conversion){
  for(q in 1:nrow(hyper.param.conversion.idx)){
    j = hyper.param.conversion.idx[q,1]
    k = hyper.param.conversion.idx[q,2]
    l = hyper.param.conversion.idx[q,3]
    conversion.ce.simulation.results.array[i,j,k,l] <- conversion.ce.simulation.results[[i]]$smoothed[q]
  }
  for(q in 1:nrow(hyper.param.conversion.ml.idx)){
    j = hyper.param.conversion.ml.idx[q,1]
    k = hyper.param.conversion.ml.idx[q,2]
    conversion.ce.simulation.ml.results.array[i,j,k] <- conversion.ce.simulation.results[[i]]$npmle[q]
  }
}

saveRDS(conversion.ce.simulation.results.array,paste0("Data/results/conversion_ce_simulations.rds"))
saveRDS(conversion.ce.simulation.ml.results.array,paste0("Data/results/conversion_ce_ml_simulations.rds"))

print("conversion cross entropy simulations complete")






## -------- Feasibility Tests -------- 


# feasibility.test.first.order.results.array <- array(NA, dim = c(n.sims.feasibility, length(h.set.feasibility), length(ker.set)))
# feasibility.test.second.order.results.array <- array(NA, dim = c(n.sims.feasibility, length(h.set.feasibility), length(ker.set)))
# 
# 
# feasibility.test.results <- simulate_experiment(RUN_PARALLEL, 
#                                                 simulation_feasibility_test_model1, 
#                                                 n.sims.feasibility)
# 
# for(i in 1:n.sims.feasibility){
#   for(q in 1:nrow(hyper.param.feasibility.idx)){
#     j = hyper.param.feasibility.idx[q,1]
#     k = hyper.param.feasibility.idx[q,2]
#     feasibility.test.first.order.results.array[i,j,k] <- feasibility.test.results[[i]][[q]]$first_order
#     feasibility.test.second.order.results.array[i,j,k] <- feasibility.test.results[[i]][[q]]$second_order
#   }
# }
# 
# saveRDS(feasibility.test.first.order.results.array,paste0("Data/results/first_order_feasibility_test_simulations.rds"))
# saveRDS(feasibility.test.second.order.results.array,paste0("Data/results/second_order_feasibility_test_simulations.rds"))
# 
# print("feasibility test simulations complete")
# 
# 
# 



## -------- Speed Tests -------- 


# speed.test.time.results.array <- array(NA, dim = c(n.sims.speed, length(h.set.speed), 2, length(mu.set)))
# speed.test.likelihood.results.array <- array(NA, dim = c(n.sims.speed, length(h.set.speed),2, length(mu.set)))
# 
# 
# speed.test.results <- simulate_experiment(RUN_PARALLEL, 
#                                           simulation_speed_test, 
#                                           n.sims.speed)
# 
# saveRDS(speed.test.results,paste0("Data/results/speed_test.rds"))
# 
#
# 
# for(i in 1:n.sims.speed){
#   for(l in 1:length(h.set.speed)){
#     for(k in 1:length(mu.set)){
#       speed.test.time.results.array[i,l,1,k] <- speed.test.results[[i]]$time[l,1,k]
#       speed.test.time.results.array[i,l,2,k] <- speed.test.results[[i]]$time[l,2,k]
#       speed.test.likelihood.results.array[i,l,1,k] <- speed.test.results[[i]]$likelihood[l,1,k]
#       speed.test.likelihood.results.array[i,l,2,k] <- speed.test.results[[i]]$likelihood[l,2,k]
#       
#     }
#   }
# }
# 
# saveRDS(speed.test.time.results.array,paste0("Data/results/speed_test_time_simulations.rds"))
# saveRDS(speed.test.likelihood.results.array,paste0("Data/results/speed_test_likelihood_simulations.rds"))
# 
# print("speed test simulations complete")
# 












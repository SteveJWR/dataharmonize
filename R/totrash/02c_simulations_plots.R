
library(gridExtra)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
#library(tidyr)
library(R.utils)


set.seed(2)

### In this R script, we produce all of the relevant plots
source("R/01_functions.R")
source("R/02a_simulations_setup.R")

##### Figure information 
png.width = 1200
png.height = 1000
png.res = 200 

title.size = 16 
axis.ticks.size = 10 
axis.size = 12  

small.title.size = 12 
small.axis.ticks.size = 8 
small.axis.size = 10 

extra.small.title.size = 8
extra.small.axis.ticks.size = 8
extra.small.axis.size = 8



### Intrinsic_variability plot model 1 

intrinsic.variability.model1.simulation.results.array <- readRDS("Data/results/intrinsic_variability_simulation_model1.rds")
iv.model1 <- summarize_simulations(intrinsic.variability.model1.simulation.results.array)

iv.model1.mean <- iv.model1$mean_results
iv.model1.se <- iv.model1$se_results

mus <- round(rep(mu.set, each = length(h.set)), 3)
hs <- rep(h.set, times = length(mu.set))

gauss.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(iv.model1.mean[,1,]), total.variation.sd = as.vector(iv.model1.se[,1,]) )
exp.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(iv.model1.mean[,2,]), total.variation.sd = as.vector(iv.model1.se[,2,]))
triangle.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(iv.model1.mean[,3,]), total.variation.sd = as.vector(iv.model1.se[,3,]))
epanechnikov.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(iv.model1.mean[,4,]), total.variation.sd = as.vector(iv.model1.se[,4,]))


# argmin is the gaussian kernel with mu = 0.01

gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
triangle.frame$mu <- factor(triangle.frame$mu)
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)







gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Gaussian Kernel") + 
  ylab("Total Variation") + 
  geom_vline(xintercept = h.model1) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))

exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Laplace Kernel") + 
  ylab("Total Variation") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))

triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Triangle Kernel") + 
  ylab("Total Variation") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))

epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Epanechnikov Kernel") + 
  ylab("Total Variation") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))



png(filename = "Plots/intrinsic_variability_simulation_model1_large.png",
    width = png.width, height = png.height, res = png.res)

ggarrange(gauss.plot, exp.plot, 
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 

png(filename = "Plots/intrinsic_variability_simulation_model1.png",
    width = 2*png.width, height = png.height, res = png.res)

ggarrange(gauss.plot, exp.plot, epanechnikov.plot,
          ncol=3, nrow=1, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 





### Intrinsic_variability plot model 2 

intrinsic.variability.model2.simulation.results.array <- readRDS("Data/results/intrinsic_variability_simulation_model2.rds")
iv.model2 <- summarize_simulations(intrinsic.variability.model2.simulation.results.array)

iv.model2.mean <- iv.model2$mean_results
iv.model2.se <- iv.model2$se_results

mus <- round(rep(mu.set, each = length(h.set)), 3)
hs <- rep(h.set, times = length(mu.set))

gauss.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(iv.model2.mean[,1,]), total.variation.sd = as.vector(iv.model2.se[,1,]) )
exp.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(iv.model2.mean[,2,]), total.variation.sd = as.vector(iv.model2.se[,2,]))
triangle.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(iv.model2.mean[,3,]), total.variation.sd = as.vector(iv.model2.se[,3,]))
epanechnikov.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(iv.model2.mean[,4,]), total.variation.sd = as.vector(iv.model2.se[,4,]))


# argmin is the gaussian kernel with mu = 0.01

gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
triangle.frame$mu <- factor(triangle.frame$mu)
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)







gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Gaussian Kernel") + 
  ylab("Total Variation") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))

exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Laplace Kernel") + 
  ylab("Total Variation") + 
  geom_vline(xintercept = h.model2) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))

triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Triangle Kernel") + 
  ylab("Total Variation") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))

epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Epanechnikov Kernel") + 
  ylab("Total Variation") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))



png(filename = "Plots/intrinsic_variability_simulation_model2_large.png",
    width = png.width, height = png.height, res = png.res)

ggarrange(gauss.plot, exp.plot, 
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 

png(filename = "Plots/intrinsic_variability_simulation_model2.png",
    width = 2*png.width, height = png.height, res = png.res)

ggarrange(gauss.plot, exp.plot, epanechnikov.plot,
          ncol=3, nrow=1, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 







#### Mu Selection -- Model 1

smoothing.selection.model1.simulation.results.array <- readRDS("Data/results/smoothing_selection_simulation_model1.rds")
# edit this to handle the case where there are only two dimensions in the array
smooth.model1 <- summarize_simulations(smoothing.selection.model1.simulation.results.array)

smooth.model1.mean <- smooth.model1$mean_results
smooth.model1.se <- smooth.model1$se_results


smoothing.model1.frame <- data.frame(log.mu = log(mu.set.long + 0.001), neg.loglike = smooth.model1.mean, neg.loglike.sd = smooth.model1.se)
smoothing.model1.plot <- ggplot(data = smoothing.model1.frame, aes(x = log.mu, y = neg.loglike)) + geom_line() + 
  geom_errorbar(aes(ymin=neg.loglike-2*neg.loglike.sd, ymax=neg.loglike+2*neg.loglike.sd), 
                width=.05, position=position_dodge(0.05)) + ylab("Negative Two Sample Log-Likelihood") + 
  ggtitle("Latent: Beta(12,5) and MKM(Gaussian,h = 2)") + xlab("log(\u00b5 + 0.001)") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  
smoothing.model1.plot



print(paste0("Optimal Mu Value: ", mu.set.long[which.min(smooth.model1.mean)]))
# 0.0193


png(filename = "Plots/smoothing_selection_simulation_model1.png",
    width = png.width, height = png.height, res = png.res)

smoothing.model1.plot

# Close the pdf file
dev.off() 





#### Mu Selection -- Model 2

smoothing.selection.model2.simulation.results.array <- readRDS("Data/results/smoothing_selection_simulation_model2.rds")
# edit this to handle the case where there are only two dimensions in the array
smooth.model2 <- summarize_simulations(smoothing.selection.model2.simulation.results.array)

smooth.model2.mean <- smooth.model2$mean_results
smooth.model2.se <- smooth.model2$se_results


smoothing.model2.frame <- data.frame(log.mu = log(mu.set.long + 0.001), neg.loglike = smooth.model2.mean, neg.loglike.sd = smooth.model2.se)
smoothing.model2.plot <- ggplot(data = smoothing.model2.frame, aes(x = log.mu, y = neg.loglike)) + geom_line() + 
  geom_errorbar(aes(ymin=neg.loglike-2*neg.loglike.sd, ymax=neg.loglike+2*neg.loglike.sd), 
                width=.05, position=position_dodge(0.05)) + ylab("Negative Two Sample Log-Likelihood") + 
  ggtitle("Latent: Beta(6,6) and MKM(Laplace,h = 1)") + xlab("log(\u00b5 + 0.001)") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  
smoothing.model2.plot



print(paste0("Optimal Mu Value: ", mu.set.long[which.min(smooth.model2.mean)]))
# 0.0139


png(filename = "Plots/smoothing_selection_simulation_model2.png",
    width = png.width, height = png.height, res = png.res)

smoothing.model2.plot

# Close the pdf file
dev.off() 






### ---- conversion cross entropy ----


conversion.ce.simulation.results.array <- readRDS("Data/results/conversion_ce_simulations.rds")
conversion.ce.simulation.ml.results.array <- readRDS("Data/results/conversion_ce_ml_simulations.rds")

conversion.ce <- summarize_simulations(conversion.ce.simulation.results.array)
conversion.ce.mean <- conversion.ce$mean_results
conversion.ce.se <- conversion.ce$se_results


conversion.ce.ml <- summarize_simulations(conversion.ce.simulation.ml.results.array)
conversion.ce.mean.ml <- conversion.ce.ml$mean_results
conversion.ce.se.ml <- conversion.ce.ml$se_results



P6 <- brewer.pal(7, "Blues")
P7 <- c(P6[2:7],"red")
mus <- round(rep(mu.set.conversion, each = length(h.set.conversion)), 4)
mus <- as.character(mus)
mus <- c(mus, rep("unregularized", length(h.set.conversion)))
hs <- rep(h.set.conversion, times = length(mu.set.conversion) + 1)

gauss.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(conversion.ce.mean[,1,]),as.vector(conversion.ce.mean.ml[,1])), cross.entropy.sd = c(as.vector(conversion.ce.se[,1,]),as.vector(conversion.ce.se.ml[,1])) )
exp.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(conversion.ce.mean[,2,]),as.vector(conversion.ce.mean.ml[,2])), cross.entropy.sd = c(as.vector(conversion.ce.se[,2,]),as.vector(conversion.ce.se.ml[,2])))



#Narrowing Colour Choices 
gauss.frame <- gauss.frame[gauss.frame$mu %in% c("0","0.0037", "0.0139", "0.0193", "0.0761","0.3","unregularized"), ]
exp.frame <- exp.frame[exp.frame$mu %in% c("0","0.0037", "0.0139", "0.0193", "0.0761","0.3","unregularized"), ]


gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)



# model selection index, exp frame 31 
gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + 
  geom_line() + #geom_errorbar(aes(ymin=cross.entropy-3*cross.entropy.sd, ymax=cross.entropy+3*cross.entropy.sd), 
                #              width=.2, position=position_dodge(0.0)) + 
  ggtitle("Gaussian Kernel") + 
  coord_cartesian(xlim = c(0,4.1), ylim = c(2.8, 3.3)) + ylab("Population Cross Entropy")  

exp.plot <- ggplot(data = exp.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + 
  geom_line() + #geom_errorbar(aes(ymin=cross.entropy-3*cross.entropy.sd, ymax=cross.entropy+3*cross.entropy.sd), 
                #              width=.2, position=position_dodge(0.0)) + 
  ggtitle("Laplace Kernel") + 
  coord_cartesian(xlim = c(0,4.1), ylim = c(2.8, 3.3)) + ylab("Population Cross Entropy")  + 
  annotate("point", x = exp.frame$h[30], y = exp.frame$cross.entropy[30], colour = "blue", shape = "x", size = 8) #+ geom_point(size = 0.6) +  geom_point(aes(x = exp.frame$h[31], y = exp.frame$cross.entropy[31], shape = "+", color = "black"))


gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + 
  geom_line() + 
  geom_errorbar(aes(ymin=cross.entropy-2*cross.entropy.sd, 
                    ymax=cross.entropy+2*cross.entropy.sd),
                width=.2, position=position_dodge(0.0)) +
  ggtitle("Gaussian Kernel")  + coord_cartesian(xlim = c(0,4.1), ylim = c(2.35, 3.0)) + 
  ylab("Population Cross Entropy") + scale_color_manual(values= as.vector(P7)) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

exp.plot <- ggplot(data = exp.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + 
  geom_line() + 
  geom_errorbar(aes(ymin=cross.entropy-2*cross.entropy.sd, 
                    ymax=cross.entropy+2*cross.entropy.sd),
                width=.2, position=position_dodge(0.0)) +
  ggtitle("Laplace Kernel") + coord_cartesian(xlim = c(0,4.1), ylim = c(2.35, 3.0)) + 
  ylab("Population Cross Entropy")  + 
  annotate("point", x = exp.frame$h[30], y = exp.frame$cross.entropy[30], colour = "blue", shape = "x", size = 8) + 
  scale_color_manual(values= as.vector(P7)) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  






png(filename = "Plots/conversion_simulation.png",
    width = 1.5*png.width, height = png.height, res = png.res)


ggarrange(gauss.plot, exp.plot,
          ncol=2, nrow=1, common.legend = TRUE, legend="right") 
# Close the pdf file
dev.off() 




### Feasibility tests 

feasibility.test.first.order.results.array <- readRDS("Data/results/first_order_feasibility_test_simulations.rds")
feasibility.test.second.order.results.array <- readRDS("Data/results/second_order_feasibility_test_simulations.rds")


gauss.frame.1 <- data.frame(matrix(NA, nrow = n.sims.feasibility, ncol = length(h.set.feasibility)))
exp.frame.1 <- data.frame(matrix(NA, nrow = n.sims.feasibility, ncol = length(h.set.feasibility)))
triangle.frame.1 <- data.frame(matrix(NA, nrow = n.sims.feasibility, ncol = length(h.set.feasibility)))
epanechnikov.frame.1 <- data.frame(matrix(NA, nrow = n.sims.feasibility, ncol = length(h.set.feasibility)))

colnames(gauss.frame.1) <- h.set.feasibility 
colnames(exp.frame.1) <- h.set.feasibility
colnames(triangle.frame.1) <- h.set.feasibility
colnames(epanechnikov.frame.1) <- h.set.feasibility


gauss.frame.2 <- data.frame(matrix(NA, nrow = n.sims.feasibility, ncol = length(h.set.feasibility)))
exp.frame.2 <- data.frame(matrix(NA, nrow = n.sims.feasibility, ncol = length(h.set.feasibility)))
triangle.frame.2 <- data.frame(matrix(NA, nrow = n.sims.feasibility, ncol = length(h.set.feasibility)))
epanechnikov.frame.2 <- data.frame(matrix(NA, nrow = n.sims.feasibility, ncol = length(h.set.feasibility)))

colnames(gauss.frame.2) <- h.set.feasibility
colnames(exp.frame.2) <- h.set.feasibility
colnames(triangle.frame.2) <- h.set.feasibility
colnames(epanechnikov.frame.2) <- h.set.feasibility




for(i in 1:n.sims.feasibility){
  gauss.frame.1[i,] <- feasibility.test.first.order.results.array[i,,1]
  exp.frame.1[i,] <- feasibility.test.first.order.results.array[i,,2]
  triangle.frame.1[i,] <- feasibility.test.first.order.results.array[i,,3]
  epanechnikov.frame.1[i,] <- feasibility.test.first.order.results.array[i,,4]
  
  gauss.frame.2[i,] <- feasibility.test.second.order.results.array[i,,1]
  exp.frame.2[i,] <- feasibility.test.second.order.results.array[i,,2]
  triangle.frame.2[i,] <- feasibility.test.second.order.results.array[i,,3]
  epanechnikov.frame.2[i,] <- feasibility.test.second.order.results.array[i,,4]
  
}


gauss.frame.1.long <- gauss.frame.1 %>% tidyr::gather(bandwidth, p.value)
exp.frame.1.long <- exp.frame.1 %>% tidyr::gather(bandwidth, p.value)
triangle.frame.1.long <- triangle.frame.1 %>% tidyr::gather(bandwidth, p.value)
epanechnikov.frame.1.long <- epanechnikov.frame.1 %>% tidyr::gather(bandwidth, p.value)

gauss.frame.1.long$bandwidth <- factor(gauss.frame.1.long$bandwidth, levels = h.set.feasibility)
exp.frame.1.long$bandwidth <- factor(exp.frame.1.long$bandwidth, levels = h.set.feasibility)
triangle.frame.1.long$bandwidth <- factor(triangle.frame.1.long$bandwidth, levels = h.set.feasibility)
epanechnikov.frame.1.long$bandwidth <- factor(epanechnikov.frame.1.long$bandwidth, levels = h.set.feasibility)

gauss.frame.2.long <- gauss.frame.2 %>% tidyr::gather(bandwidth, p.value)
exp.frame.2.long <- exp.frame.2 %>% tidyr::gather(bandwidth, p.value)
triangle.frame.2.long <- triangle.frame.2 %>% tidyr::gather(bandwidth, p.value)
epanechnikov.frame.2.long <- epanechnikov.frame.2 %>% tidyr::gather(bandwidth, p.value)

gauss.frame.2.long$bandwidth <- factor(gauss.frame.2.long$bandwidth, levels = h.set.feasibility)
exp.frame.2.long$bandwidth <- factor(exp.frame.2.long$bandwidth, levels = h.set.feasibility)
triangle.frame.2.long$bandwidth <- factor(triangle.frame.2.long$bandwidth, levels = h.set.feasibility)
epanechnikov.frame.2.long$bandwidth <- factor(epanechnikov.frame.2.long$bandwidth, levels = h.set.feasibility)


p.gauss.1 <- ggplot(data = gauss.frame.1.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Gaussian Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  

p.exp.1 <- ggplot(data = exp.frame.1.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Laplace Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  

p.triangle.1 <- ggplot(data = triangle.frame.1.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Triangle Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  

p.epanechnikov.1 <- ggplot(data = epanechnikov.frame.1.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Epanechnikov Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  



p.gauss.2 <- ggplot(data = gauss.frame.2.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Gaussian Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  

p.exp.2 <- ggplot(data = exp.frame.2.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Laplace Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  

p.triangle.2 <- ggplot(data = triangle.frame.2.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Triangle Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  

p.epanechnikov.2 <- ggplot(data = epanechnikov.frame.2.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Epanechnikov Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  



grid.arrange(arrangeGrob(p.gauss.1, p.exp.1, ncol = 2),                             
             arrangeGrob(p.triangle.1, p.epanechnikov.1, ncol = 2), 
             nrow = 2)   

grid.arrange(arrangeGrob(p.gauss.2, p.exp.2, ncol = 2),                             
             arrangeGrob(p.triangle.2, p.epanechnikov.2, ncol = 2), 
             nrow = 2)   



png(filename = "Plots/first_order_feasibililty_simulations_large.png",
    width = png.width, height = png.height, res = png.res)

ggarrange(p.gauss.1, p.exp.1, p.triangle.1, 
          p.epanechnikov.1, nrow = 2, ncol = 2, 
          common.legend = TRUE, legend="right")   

# Close the pdf file
dev.off() 


png(filename = "Plots/first_order_feasibililty_simulations.png",
    width = 2*png.width, height = png.height, res = png.res)

ggarrange(p.gauss.1, p.exp.1, p.epanechnikov.1, nrow = 1, ncol = 3, 
          common.legend = TRUE, legend="right")   

# Close the pdf file
dev.off() 

png(filename = "Plots/second_order_feasibililty_simulations_large.png",
    width = png.width, height = png.height, res = png.res)

ggarrange(p.gauss.2, p.exp.2, p.triangle.2, 
          p.epanechnikov.2, nrow = 2, ncol = 2, 
          common.legend = TRUE, legend="right")   

# Close the pdf file
dev.off() 

png(filename = "Plots/second_order_feasibililty_simulations.png",
    width = 2*png.width, height = png.height, res = png.res)

ggarrange(p.gauss.2, p.exp.2, p.epanechnikov.2, nrow = 1, ncol = 3, 
          common.legend = TRUE, legend="right")   

# Close the pdf file
dev.off() 




### Speed tests 

speed.test.time.results.array <- readRDS("Data/results/speed_test_time_simulations.rds")
speed.test.likelihood.results.array <- readRDS("Data/results/speed_test_likelihood_simulations.rds")
speed.test.likelihood.rel.diff.results.array <- (speed.test.likelihood.results.array[, ,1, ] - speed.test.likelihood.results.array[, ,2, ])/(abs(speed.test.likelihood.results.array[, ,1, ])) 

speed.test.time.summary <- summarize_simulations(speed.test.time.results.array)
speed.test.time.mean <- speed.test.time.summary$mean_results
speed.test.time.se <- speed.test.time.summary$se_results

speed.test.like.summary <- summarize_simulations(speed.test.likelihood.rel.diff.results.array)
speed.test.like.mean <- speed.test.like.summary$mean_results
speed.test.like.se <- speed.test.like.summary$se_results



mus <- round(rep(mu.set, each = length(h.set.speed)), 3)
hs <- rep(h.set.speed, times = length(mu.set))

npem.speed.frame <- data.frame(h = hs, mu = mus, time = as.vector(speed.test.time.mean[,1,]), sd = as.vector(speed.test.time.se[,1,]) )
gp.speed.frame <- data.frame(h = hs, mu = mus, time = as.vector(speed.test.time.mean[,2,]), sd = as.vector(speed.test.time.se[,2,]) )
diff.likelihood.frame <- data.frame(h = hs, mu = mus, likelihood = as.vector(speed.test.like.mean[,]), sd = as.vector(speed.test.like.se[,]))


npem.speed.frame$mu <- factor(npem.speed.frame$mu)
gp.speed.frame$mu <- factor(gp.speed.frame$mu)
diff.likelihood.frame$mu <- factor(diff.likelihood.frame$mu)



y.min.like = min(c(speed.test.like.mean[,]))
y.max.like = max(c(speed.test.like.mean[,]))

npem.speed.plot <- ggplot(data = npem.speed.frame, aes(x = log(h), y = time, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=time-2*sd, ymax=time+2*sd), width=.2,
                              position=position_dodge(0.05)) + ggtitle("NPEM Speed") + ylab("Fit Time (s)") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  

gp.speed.plot <- ggplot(data = gp.speed.frame, aes(x = log(h), y = time, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=time-2*sd, ymax=time+2*sd), width=.2,
                              position=position_dodge(0.05)) + ggtitle("GP Speed") + ylab("Fit Time (s)") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  


diff.likelihood.plot <- ggplot(data = diff.likelihood.frame, aes(x = log(h), y = likelihood, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=likelihood-2*sd, ymax=likelihood+2*sd), width=.2,
                              position=position_dodge(0.05)) + ggtitle("Log-Likelihood Difference (NPEM - GP)/|NPEM|") + 
  theme(axis.text.x = element_text( size = small.axis.ticks.size),
        axis.text.y = element_text( size = small.axis.ticks.size),  
        axis.title.x = element_text( size = small.axis.size),
        axis.title.y = element_text( size = small.axis.size),
        title = element_text( size = small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=small.axis.size), #change legend title font size
        legend.text = element_text(size=small.axis.ticks.size))  




grid.arrange(arrangeGrob(npem.speed.plot, gp.speed.plot, ncol = 2),                             
             diff.likelihood.plot, 
             nrow = 2)   

png(filename = "plots/speed_test_simulations.png",
    width = png.width, height = png.height, res = png.res)


grid.arrange(arrangeGrob(npem.speed.plot, gp.speed.plot, ncol = 2),                             
             diff.likelihood.plot, 
             nrow = 2)   


# Close the pdf file
dev.off() 




#### Lastly we produce an example of the effect of regularization on a mixing density problem


sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                          beta1.y = beta1.model1, beta2.y = beta2.model1, 
                          cond.y = cond.model1, pair.obs = F)
train.p.hat <- compute_edf(sim.data[,1], N)

sim.data.large <- simulate_beta(n.ind = 100000, n.obs.per.ind = 2, 
                                beta1.y = beta1.model1, beta2.y = beta2.model1, 
                                cond.y = cond.model1, pair.obs = F)
true.p.hat <- compute_edf(sim.data.large[,1], N)





model.estimate0 <- estimate_mixing_numeric(p.hat = train.p.hat, 
                                           A.matrix = A.matrix.model1, 
                                           mu = 0)

model.estimate0.001 <- estimate_mixing_numeric(p.hat = train.p.hat, 
                                          A.matrix = A.matrix.model1, 
                                          mu = 0.001)

model.estimate0.01 <- estimate_mixing_numeric(p.hat = train.p.hat, 
                                         A.matrix = A.matrix.model1, 
                                         mu = 0.01)

model.estimate0.1 <- estimate_mixing_numeric(p.hat = train.p.hat, 
                                        A.matrix = A.matrix.model1, 
                                        mu = 0.1) 

model.estimate0.4 <- estimate_mixing_numeric(p.hat = train.p.hat, 
                                        A.matrix = A.matrix.model1, 
                                        mu = 0.4) 

gamma.seq <- seq(0,1,length.out = R_bins)

cdf.0 <- rep(NA, R_bins)
cdf.0.001 <- rep(NA, R_bins)
cdf.0.01 <- rep(NA, R_bins)
cdf.0.1 <- rep(NA, R_bins)
cdf.0.4 <- rep(NA, R_bins)

for(k in 1:R_bins){
  cdf.0[k] <- sum(model.estimate0$latent[1:k])
  cdf.0.001[k] <- sum(model.estimate0.001$latent[1:k])
  cdf.0.01[k] <- sum(model.estimate0.01$latent[1:k])
  cdf.0.1[k] <- sum(model.estimate0.1$latent[1:k])
  cdf.0.4[k] <- sum(model.estimate0.4$latent[1:k])
}


latent.true <- dbeta(gamma.seq, beta1.model1, beta2.model1)
cdf.true <- pbeta(gamma.seq, beta1.model1, beta2.model1)

y.seq = 0:N



### Observed Data
obs.mu0 <- data.frame(y = rep(y.seq, 3), 
                      labels = rep(c("True", "Empirical", "Estimate"), 
                                   each = length(y.seq)),
                      Probability = c(true.p.hat, 
                                      train.p.hat, 
                                      as.vector(model.estimate0$observed)))


obs.mu0.001 <- data.frame(y = rep(y.seq, 3), 
                          labels = rep(c("True", "Empirical", "Estimate"), 
                                       each = length(y.seq)),
                          Probability = c(true.p.hat, 
                                          train.p.hat, 
                                          as.vector(model.estimate0.001$observed)))


obs.mu0.01 <- data.frame(y = rep(y.seq, 3), 
                         labels = rep(c("True", "Empirical", "Estimate"), 
                                      each = length(y.seq)),
                         Probability = c(true.p.hat, 
                                         train.p.hat, 
                                         as.vector(model.estimate0.01$observed)))


obs.mu0.1 <- data.frame(y = rep(y.seq, 3), 
                        labels = rep(c("True", "Empirical", "Estimate"), 
                                     each = length(y.seq)),
                        Probability = c(true.p.hat, 
                                        train.p.hat, 
                                        as.vector(model.estimate0.1$observed)))

obs.mu0.4 <- data.frame(y = rep(y.seq, 3), 
                        labels = rep(c("True", "Empirical", "Estimate"), 
                                     each = length(y.seq)),
                        Probability = c(true.p.hat, 
                                        train.p.hat, 
                                        as.vector(model.estimate0.4$observed)))


### Latent Density 


dens.mu0 <- data.frame(y = rep(gamma.seq, 2), 
                       labels = rep(c("True", "Estimate"), 
                                    each = length(gamma.seq)),
                       density = c(latent.true, 
                                   R_bins*as.vector(model.estimate0$latent)))


dens.mu0.001 <- data.frame(y = rep(gamma.seq, 2), 
                           labels = rep(c("True", "Estimate"), 
                                        each = length(gamma.seq)),
                           density = c(latent.true, 
                                       R_bins*as.vector(model.estimate0.001$latent)))

dens.mu0.01 <- data.frame(y = rep(gamma.seq, 2), 
                          labels = rep(c("True", "Estimate"), 
                                       each = length(gamma.seq)),
                          density = c(latent.true, 
                                      R_bins*as.vector(model.estimate0.01$latent)))

dens.mu0.1 <- data.frame(y = rep(gamma.seq, 2), 
                         labels = rep(c("True", "Estimate"), 
                                      each = length(gamma.seq)),
                         density = c(latent.true, 
                                     R_bins*as.vector(model.estimate0.1$latent)))

dens.mu0.4 <- data.frame(y = rep(gamma.seq, 2), 
                         labels = rep(c("True", "Estimate"), 
                                      each = length(gamma.seq)),
                         density = c(latent.true, 
                                     R_bins*as.vector(model.estimate0.4$latent)))




### CDF 


cdf.mu0 <-  data.frame(y = rep(gamma.seq, 2), 
                       labels = rep(c("True", "Estimate"), 
                                    each = length(gamma.seq)),
                       cdf = c(cdf.true, 
                               cdf.0))

cdf.mu0.001 <- data.frame(y = rep(gamma.seq, 2), 
                          labels = rep(c("True", "Estimate"), 
                                       each = length(gamma.seq)),
                          cdf = c(cdf.true, 
                                  cdf.0.001))

cdf.mu0.01 <- data.frame(y = rep(gamma.seq, 2), 
                         labels = rep(c("True", "Estimate"), 
                                      each = length(gamma.seq)),
                         cdf = c(cdf.true,
                                 cdf.0.01))

cdf.mu0.1 <- data.frame(y = rep(gamma.seq, 2), 
                        labels = rep(c("True", "Estimate"), 
                                     each = length(gamma.seq)),
                        cdf = c(cdf.true,
                                cdf.0.1))

cdf.mu0.4 <- data.frame(y = rep(gamma.seq, 2), 
                        labels = rep(c("True", "Estimate"), 
                                     each = length(gamma.seq)),
                        cdf = c(cdf.true, 
                                cdf.0.4))




p.obs.0 <- ggplot(data = obs.mu0, aes(x = y, 
                                      y = Probability, 
                                      group = labels, 
                                      color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0)) + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) 


p.obs.0.001 <- ggplot(data = obs.mu0.001, aes(x = y, 
                                              y = Probability, 
                                              group = labels, 
                                              color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0.001)) + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) 



p.obs.0.01 <- ggplot(data = obs.mu0.01, aes(x = y, 
                                            y = Probability, 
                                            group = labels, 
                                            color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0.01)) + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) 


p.obs.0.1 <- ggplot(data = obs.mu0.1, aes(x = y, 
                                          y = Probability, 
                                          group = labels, 
                                          color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0.1)) + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) 


p.obs.0.4 <- ggplot(data = obs.mu0.4, aes(x = y, 
                                          y = Probability, 
                                          group = labels, 
                                          color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0.4)) + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) 



##### Latent


p.dens.0 <- ggplot(data = dens.mu0, aes(x = y, 
                                        y = density, 
                                        group = labels, 
                                        color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) #+ 
# labs(col="Density") + annotate(geom="text", x=0.1, y= 500, label=paste0("\u00b5 = ", 0),
#                                color="blue", size = 6)


p.dens.0


p.dens.0.001 <- ggplot(data = dens.mu0.001, aes(x = y, 
                                                y = density, 
                                                group = labels, 
                                                color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0.001)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) #+ 
# labs(col="Density") + annotate(geom="text", x=0.1, y= 45, label=paste0("\u00b5 = ", 0.001),
#                                color="blue", size = 6)

p.dens.0.001



p.dens.0.01 <- ggplot(data = dens.mu0.01, aes(x = y, 
                                              y = density, 
                                              group = labels, 
                                              color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0.01)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) #+ 
# labs(col="Density") + annotate(geom="text", x=0.1, y= 12, label=paste0("\u00b5 = ", 0.01),
#                                color="blue", size = 6)

p.dens.0.01


p.dens.0.1 <- ggplot(data = dens.mu0.1, aes(x = y, 
                                            y = density, 
                                            group = labels, 
                                            color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0.1)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) #+ 
# labs(col="Density") + annotate(geom="text", x=0.1, y= 6, label=paste0("\u00b5 = ", 0.1),
#                                color="blue", size = 6)

p.dens.0.1

p.dens.0.4 <- ggplot(data = dens.mu0.4, aes(x = y, 
                                            y = density, 
                                            group = labels, 
                                            color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0.4)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) #+ 
# labs(col="Density") + annotate(geom="text", x=0.1, y= 4, label=paste0("\u00b5 = ", 0.4),
#            color="blue", size = 6)

p.dens.0.4


##### CDF


p.cdf.0 <- ggplot(data = cdf.mu0, aes(x = y, 
                                      y = cdf, 
                                      group = labels, 
                                      color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) 


p.cdf.0.001 <- ggplot(data = cdf.mu0.001, aes(x = y, 
                                              y = cdf, 
                                              group = labels, 
                                              color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0.001)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) 



p.cdf.0.01 <- ggplot(data = cdf.mu0.01, aes(x = y, 
                                            y = cdf, 
                                            group = labels, 
                                            color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0.01)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) 


p.cdf.0.1 <- ggplot(data = cdf.mu0.1, aes(x = y, 
                                          y = cdf, 
                                          group = labels, 
                                          color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0.1)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size)) 


p.cdf.0.4 <- ggplot(data = cdf.mu0.4, aes(x = y, 
                                          y = cdf, 
                                          group = labels, 
                                          color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0.4)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = extra.small.axis.ticks.size),
        axis.text.y = element_text( size = extra.small.axis.ticks.size),  
        axis.title.x = element_text( size = extra.small.axis.size),
        axis.title.y = element_text( size = extra.small.axis.size),
        title = element_text( size = extra.small.title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=extra.small.axis.size), #change legend title font size
        legend.text = element_text(size=extra.small.axis.ticks.size))  






png(filename = "Plots/Example Regularization.png",
    width = (1.5)*png.width, height = png.height, res = png.res)


ggarrange(p.obs.0, p.dens.0, p.cdf.0,
          p.obs.0.001, p.dens.0.001, p.cdf.0.001,
          p.obs.0.01, p.dens.0.01, p.cdf.0.01,
          p.obs.0.1, p.dens.0.1, p.cdf.0.1,
          p.obs.0.4, p.dens.0.4, p.cdf.0.4,
          ncol=3, nrow=5, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 









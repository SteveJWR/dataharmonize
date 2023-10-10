
source("R/01_functions.R")
source("R/03a_application_setup.R")

library(RColorBrewer)
library(colorRamps)
library(pals)

cond.y.model.0 <- readRDS("Data/results/conditional_model_selection_mu0_mmse.rds")
cond.y.model.0.001 <- readRDS("Data/results/conditional_model_selection_mu0001_mmse.rds")
cond.y.model.0.1 <- readRDS("Data/results/conditional_model_selection_mu01_mmse.rds")

cond.z.model.0 <- readRDS("Data/results/conditional_model_selection_mu0_moca.rds")
cond.z.model.0.001 <- readRDS("Data/results/conditional_model_selection_mu0001_moca.rds")
cond.z.model.0.1 <- readRDS("Data/results/conditional_model_selection_mu01_moca.rds")


y.smooth <- readRDS("Data/results/smoothing_selection_mmse.rds")
z.smooth <- readRDS("Data/results/smoothing_selection_moca.rds")

conversion.results.array <- readRDS("Data/results/regularized_cross_entropy.rds")
conversion.results.array.ml <- readRDS("Data/results/unregularized_cross_entropy.rds")
conversion.results.array.para <- readRDS("Data/results/parametric_cross_entropy.rds")
ce.zscore <- readRDS("Data/results/naive_cross_entropy.rds")




##### Figure information 
png.width = 1200
png.height = 900
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






#### Intrinsic Variability MMSE Plots



mus <- round(rep(mu.set, each = length(h.set)), 3)
hs <- rep(h.set, times = length(ker.set))

gauss.tv <- c(cond.y.model.0$results[1:length(h.set)],
              cond.y.model.0.001$results[1:length(h.set)],
              cond.y.model.0.1$results[1:length(h.set)])

exp.tv <- c(cond.y.model.0$results[(length(h.set) + 1:length(h.set))],
            cond.y.model.0.001$results[(length(h.set) + 1:length(h.set))],
            cond.y.model.0.1$results[(length(h.set) + 1:length(h.set))])

triangle.tv <- c(cond.y.model.0$results[(2*length(h.set) + 1:length(h.set))],
                 cond.y.model.0.001$results[(2*length(h.set) + 1:length(h.set))],
                 cond.y.model.0.1$results[(2*length(h.set) + 1:length(h.set))])

epanechnikov.tv <- c(cond.y.model.0$results[(3*length(h.set) + 1:length(h.set))],
                     cond.y.model.0.001$results[(3*length(h.set) + 1:length(h.set))],
                     cond.y.model.0.1$results[(3*length(h.set) + 1:length(h.set))])


gauss.frame <- data.frame(h = h.set, mu = mus, total.variation = gauss.tv)
exp.frame <- data.frame(h = h.set, mu = mus, total.variation = exp.tv)
triangle.frame <- data.frame(h = h.set, mu = mus, total.variation = triangle.tv)
epanechnikov.frame <- data.frame(h = h.set, mu = mus, total.variation = epanechnikov.tv)


binom.tv <- cond.y.model.0.001$results[(4*length(h.set) + 1)]

# argmin is the gaussian kernel with mu = 0.01

gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
triangle.frame$mu <- factor(triangle.frame$mu)
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)







gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Gaussian Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                   axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                   axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                   axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                   title = element_text( size = title.size),
                                                                                                                                                                                                                   legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                   legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                   legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                   legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                   legend.text = element_text(size=axis.ticks.size))
exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Laplace Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                              axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                              axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                              axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                              title = element_text( size = title.size),
                                                                                                                                                                                                              legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                              legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                              legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                              legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                              legend.text = element_text(size=axis.ticks.size))
triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Triangle Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                          axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                          axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                          axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                          title = element_text( size = title.size),
                                                                                                                                                                                                                          legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                          legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                          legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                          legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                          legend.text = element_text(size=axis.ticks.size))
epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Epanechnikov Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                      axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                      axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                      axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                      title = element_text( size = title.size),
                                                                                                                                                                                                                                      legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                      legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                      legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                      legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                      legend.text = element_text(size=axis.ticks.size))



png(filename = "Plots/Intrinsic_Variability_MMSE.png",
    width = png.width, height = png.height, res = png.res)

ggarrange(gauss.plot, exp.plot,
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 



#### Intrinsic Variability MOCA Plots

mus <- round(rep(mu.set, each = length(h.set)), 3)
hs <- rep(h.set, times = length(ker.set))

gauss.tv <- c(cond.z.model.0$results[1:length(h.set)],
              cond.z.model.0.001$results[1:length(h.set)],
              cond.z.model.0.1$results[1:length(h.set)])

exp.tv <- c(cond.z.model.0$results[(length(h.set) + 1:length(h.set))],
            cond.z.model.0.001$results[(length(h.set) + 1:length(h.set))],
            cond.z.model.0.1$results[(length(h.set) + 1:length(h.set))])

triangle.tv <- c(cond.z.model.0$results[(2*length(h.set) + 1:length(h.set))],
                 cond.z.model.0.001$results[(2*length(h.set) + 1:length(h.set))],
                 cond.z.model.0.1$results[(2*length(h.set) + 1:length(h.set))])

epanechnikov.tv <- c(cond.z.model.0$results[(3*length(h.set) + 1:length(h.set))],
                     cond.z.model.0.001$results[(3*length(h.set) + 1:length(h.set))],
                     cond.z.model.0.1$results[(3*length(h.set) + 1:length(h.set))])


gauss.frame <- data.frame(h = h.set, mu = mus, total.variation = gauss.tv)
exp.frame <- data.frame(h = h.set, mu = mus, total.variation = exp.tv)
triangle.frame <- data.frame(h = h.set, mu = mus, total.variation = triangle.tv)
epanechnikov.frame <- data.frame(h = h.set, mu = mus, total.variation = epanechnikov.tv)


binom.tv <- cond.z.model.0.001$results[(4*length(h.set) + 1)]

# argmin is the gaussian kernel with mu = 0.01

gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
triangle.frame$mu <- factor(triangle.frame$mu)
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)







gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Gaussian Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                   axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                   axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                   axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                   title = element_text( size = title.size),
                                                                                                                                                                                                                   legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                   legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                   legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                   legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                   legend.text = element_text(size=axis.ticks.size))
exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Laplace Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                              axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                              axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                              axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                              title = element_text( size = title.size),
                                                                                                                                                                                                              legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                              legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                              legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                              legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                              legend.text = element_text(size=axis.ticks.size))
triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Triangle Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                          axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                          axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                          axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                          title = element_text( size = title.size),
                                                                                                                                                                                                                          legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                          legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                          legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                          legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                          legend.text = element_text(size=axis.ticks.size))
epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Epanechnikov Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                      axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                      axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                      axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                      title = element_text( size = title.size),
                                                                                                                                                                                                                                      legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                      legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                      legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                      legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                      legend.text = element_text(size=axis.ticks.size))



png(filename = "Plots/Intrinsic_Variability_MOCA.png",
    width = png.width, height = png.height, res = png.res)

ggarrange(gauss.plot, exp.plot,
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 







#### Smoothing parameter selection: 


mu.frame <- data.frame(log.mu = log(mu.set.smooth + 0.001), neg.likelihood = -y.smooth$results)
mu.plot <- ggplot(data = mu.frame, aes(x = log.mu, y = neg.likelihood)) + geom_line() + ylab("Negative Two Sample Log-Likelihood")+ ggtitle("Regularization With Selected Laplace Model") + xlab("log(\u00b5 + 0.001)") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                title = element_text( size = title.size),
                                                                                                                                                                                                                                legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                legend.text = element_text(size=axis.ticks.size))
mu.plot


png(filename = "Plots/smoothing_selection_MMSE.png",
    width = png.width, height = png.height, res = png.res)

mu.plot
# Close the pdf file
dev.off() 






mu.frame <- data.frame(log.mu = log(mu.set.smooth + 0.001), neg.likelihood = -z.smooth$results)
mu.plot <- ggplot(data = mu.frame, aes(x = log.mu, y = neg.likelihood)) + geom_line() + ylab("Negative Two Sample Log-Likelihood")+ ggtitle("Regularization With Selected Binomial Model") + xlab("log(\u00b5 + 0.001)") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                title = element_text( size = title.size),
                                                                                                                                                                                                                                legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                legend.text = element_text(size=axis.ticks.size))
mu.plot


png(filename = "Plots/smoothing_selection_MOCA.png",
    width = png.width, height = png.height, res = png.res)

mu.plot
# Close the pdf file
dev.off() 




#### Plotting the conversion cross entropy. 




mus <- round(rep(mu.set.conversion, each = length(h.set)), 5)
mus <- as.character(mus)
mus <- c(mus, rep("Binomial", length(h.set)), rep("Parametric", length(h.set)), rep("Unregularized", length(h.set)), rep("Z-Score Matching", length(h.set)))
hs <- rep(h.set, times = length(mu.set.conversion) + 4)

binom.results <- conversion.results.array[61,]

gauss.idx <- 1:length(h.set)
exp.idx <- 1:length(h.set) + length(h.set)
triangle.idx <- 1:length(h.set) + 2*length(h.set)
epanechnikov.idx <- 1:length(h.set) + 3*length(h.set)



gauss.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(conversion.results.array[gauss.idx,]), rep(min(binom.results), length(h.set)),as.vector(conversion.results.array.para[gauss.idx]) ,as.vector(conversion.results.array.ml[gauss.idx]), rep(ce.zscore, length(h.set)) ) )
exp.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(conversion.results.array[exp.idx,]), rep(min(binom.results), length(h.set)),as.vector(conversion.results.array.para[exp.idx]) ,as.vector(conversion.results.array.ml[exp.idx]), rep(ce.zscore, length(h.set)) ) )
triangle.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(conversion.results.array[triangle.idx,]), rep(min(binom.results), length(h.set)),as.vector(conversion.results.array.para[triangle.idx]) ,as.vector(conversion.results.array.ml[triangle.idx]), rep(ce.zscore, length(h.set)) ) )
epanechnikov.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(conversion.results.array[epanechnikov.idx,]), rep(min(binom.results), length(h.set)),as.vector(conversion.results.array.para[epanechnikov.idx]) ,as.vector(conversion.results.array.ml[epanechnikov.idx]), rep(ce.zscore, length(h.set)) ) )




gauss.frame$mu <- factor(gauss.frame$mu, levels = as.character(c(round(mu.set.conversion, 5), "Binomial", "Parametric", "Unregularized", "Z-Score Matching")))
exp.frame$mu <- factor(exp.frame$mu, levels = as.character(c(round(mu.set.conversion, 5), "Binomial", "Parametric", "Unregularized", "Z-Score Matching")))
triangle.frame$mu <- factor(triangle.frame$mu, levels = as.character(c(round(mu.set.conversion, 5), "Binomial", "Parametric", "Unregularized", "Z-Score Matching")))
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu, levels = as.character(c(round(mu.set.conversion, 5), "Binomial", "Parametric", "Unregularized", "Z-Score Matching")))


mus.small <- round(rep(mu.set.conversion, each = length(h.set)), 5)
mus.small <- as.character(mus.small)
mus.small <- c(mus.small, rep("Parametric", length(h.set)), rep("Unregularized", length(h.set)), rep("Z-Score Matching", length(h.set)))
hs.small <- rep(h.set, times = length(mu.set.conversion) + 3)
exp.frame.small <- exp.frame




exp.frame.small <- exp.frame.small[exp.frame.small$mu %in% c("0","0.00043","0.00127","0.01129","0.0336","0.1","Parametric","Unregularized","Z-Score Matching"), ]


exp.frame.small$mu <- factor(exp.frame.small$mu, levels = c("0","0.00043","0.00127","0.01129","0.0336","0.1","Parametric","Unregularized","Z-Score Matching"))


gauss.frame.small <- gauss.frame
gauss.frame.small <- gauss.frame.small[gauss.frame.small$mu %in% c("0","0.00043","0.00127","0.01129","0.0336","0.1","Parametric","Unregularized","Z-Score Matching"), ]


gauss.frame.small$mu <- factor(gauss.frame.small$mu, levels = c("0","0.00043","0.00127","0.01129","0.0336","0.1","Parametric","Unregularized","Z-Score Matching"))






P6 <- brewer.pal(7, "Blues")
P6 <- P6[2:7]
P9 <- c(P6,"green","red","orange")

P12 <- c(brewer.pal(9, "Blues")[2:9],"green","red","orange","purple")
gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + 
  geom_line()  + ggtitle("Gaussian Kernel")  + 
  ylab("Crosswalk Sample Cross Entropy") + 
  coord_cartesian(xlim = c(0,4), ylim = c(950, 1100)) + 
  scale_color_manual(values= as.vector(P12)) + 
  annotate("point", x = gauss.frame$h[82], y = gauss.frame$cross.entropy[82], colour = "red", shape = "x", size = 6) + 
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
  geom_line()  + ggtitle("Laplace Kernel")  + ylab("Crosswalk Sample Cross Entropy")  + 
  coord_cartesian(xlim = c(0,4), ylim = c(950, 1100))  + 
  scale_color_manual(values= as.vector(P12)) + 
  annotate("point", x = exp.frame$h[80], y = exp.frame$cross.entropy[80], colour = "blue", shape = "x", size = 6) + 
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

triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + 
  geom_line()  + ggtitle("Triangle Kernel")  + 
  ylab("Crosswalk Sample Cross Entropy")  + 
  coord_cartesian(xlim = c(0,12), ylim = c(950, 1100)) + 
  scale_color_manual(values= as.vector(P12)) + 
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

epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + 
  geom_line()  + ggtitle("Epanechnikov Kernel")  + 
  ylab("Crosswalk Sample Cross Entropy")  + 
  coord_cartesian(xlim = c(0,12), ylim = c(950, 1100)) + 
  scale_color_manual(values= as.vector(P12)) + 
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
# argmin is the gaussian kernel with mu = 0.01






gauss.plot.small <- ggplot(data = gauss.frame.small, aes(x = h, y = cross.entropy, color = mu, group = mu)) + 
  geom_line()  + ggtitle("Gaussian Kernel") + 
  ylab("Crosswalk Sample Cross Entropy") + 
  coord_cartesian(xlim = c(0,4), ylim = c(950, 1100)) + 
  scale_color_manual(values= as.vector(P9)) + 
  annotate("point", x = gauss.frame.small$h[52], y = gauss.frame.small$cross.entropy[52], colour = "red", shape = "x", size = 8) + 
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




gauss.plot.small




exp.plot.small <- ggplot(data = exp.frame.small, aes(x = h, y = cross.entropy, color = mu, group = mu)) + 
  geom_line()  + ggtitle("Laplace Kernel") + 
  ylab("Crosswalk Sample Cross Entropy") + 
  coord_cartesian(xlim = c(0,4), ylim = c(950, 1100)) + 
  scale_color_manual(values= as.vector(P9)) + 
  annotate("point", x = exp.frame.small$h[50], y = exp.frame.small$cross.entropy[50], colour = "blue", shape = "x", size = 8) + 
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
exp.plot.small



gauss.plot  
exp.plot 
triangle.plot 
epanechnikov.plot

grid.arrange(arrangeGrob(gauss.plot, exp.plot, ncol = 2),
             arrangeGrob(triangle.plot, epanechnikov.plot, ncol = 2),
             nrow = 2)   


# Extra Wide Plot
png(filename = "Plots/conversion_NACC_MMSE_MOCA.png",
    width = (1.5)*png.width, height = png.height, res = png.res)

ggarrange(gauss.plot.small, exp.plot.small,
          ncol=2, nrow=1, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 


png(filename = "Plots/conversion_NACC_MMSE_MOCA_large.png",
    width = 2*png.width, height = 2*png.height, res = png.res)


ggarrange(gauss.plot, exp.plot,
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 






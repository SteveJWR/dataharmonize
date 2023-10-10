




X.ref = sim.data[idx1,]
y.ref = sim.data$Y[idx1]
Y.train
Z.train
cond.y
cond.z

mu.y/n
mu.z/n
ref.cols
ker.set
R.bins = 1000
threshold = 5 *
  10^(-5)
max.iter = 50
verbose = F





#uses bivariate data






#TODO: Check if first and second order feasibility tests are included. tryt



# saveRDS(mixture.y, "data/mixture_y.rds")
# saveRDS(mixture.z, "data/mixture_z.rds")
mixture.y <- readRDS("data/mixture_y.rds")
mixture.z <- readRDS("data/mixture_z.rds")



library(reshape2)
library(ggplot2)
library(Matrix)


p.gz <- as.matrix(p.y.gamma.joint2)
longData<-melt(as.matrix(p.y.gamma.joint2))
#longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))


plot(rowMeans(p.gz.coarse - p.gz.fine))

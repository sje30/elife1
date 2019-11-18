# Figure 3: Promoter strength, YFP, fitness

promoter_yfp_data <- read_excel("promoter-yfp data.xlsx", 
                                col_types = c("text", "numeric", "numeric", 
                                              "numeric"))


require(ggplot2)
require(data.table)
library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

prom.data <- promoter_yfp_data
prom.data <- prom.data[seq(dim(prom.data)[1],1),]
long_format <- melt(prom.data, id = "Mutation")
levels(long_format$variable) <- c("Change in TATA box strength",
                                  "YFP expression",
                                  "Relative fitness")


p <- ggplot(long_format, aes(x = Mutation,
                             y = value-1,
                             fill = variable))


p + geom_bar(stat = "summary",
             fun.y = "mean",
             width = 0.8,
             position = position_dodge(width = 0.8),
             colour = "black",
             size = 0.2) + 
  
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.8),
               width = 0.5) +
  geom_hline(yintercept = 0) +
  xlab("Mutation\n") + ylab("\nRelative units") +
  guides(fill = guide_legend(title = "Variable", reverse = TRUE)) +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.35), 
        legend.text = element_text(size = 10),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11)) +
  coord_flip()




### Promoter strength analysis
# use permutation of residuals in groups

## remove the data with missing promoter strength estimates
promoter_yfp_data = promoter_yfp_data[!is.na(promoter_yfp_data$Delta.promoter.strength),]


### YFP expression vs promoter strength: can promoter strength predict YFP fluorescence?
N=10000
T.dist=NULL
for (n in 1:N){
  new.order = sample(x = unique(promoter_yfp_data$Delta.promoter.strength), size = length(unique(promoter_yfp_data$Delta.promoter.strength)), replace = F)
  sample.Delta.promoter.strength = rep(new.order, each = 6)
  sample.mod = summary(lm(YFP.expression ~ sample.Delta.promoter.strength, data = promoter_yfp_data))
  T.dist[n] = sample.mod$coefficients[6]
}
real.T = summary(lm(YFP.expression ~ Delta.promoter.strength, data = promoter_yfp_data))$coefficients[6]
(P.value = sum(real.T<T.dist)/N)
# P = 0.1668

plot(promoter_yfp_data$Delta.promoter.strength, promoter_yfp_data$YFP.expression)
summary(mod<-lm(YFP.expression~Delta.promoter.strength, data = promoter_yfp_data))
abline(mod)

# polynomial fit and significance
summary(mod1<-lm(YFP.expression~Delta.promoter.strength, data = mean.data))
summary(mod2<-lm(YFP.expression~Delta.promoter.strength+I(Delta.promoter.strength^2), data = mean.data))
lrtest(mod1, mod2)

### fitness vs promoter strength
N=10000
T.dist=NULL
for (n in 1:N){
  new.order = sample(x = unique(promoter_yfp_data$Delta.promoter.strength), size = length(unique(promoter_yfp_data$Delta.promoter.strength)), replace = F)
  sample.Delta.promoter.strength = rep(new.order, each = 6)
  sample.mod = summary(lm(Fitness ~ sample.Delta.promoter.strength, data = promoter_yfp_data))
  # sample.mod = summary(lm(Fitness ~ log(sample.Delta.promoter.strength), data = promoter_yfp_data))
  T.dist[n] = sample.mod$coefficients[6]
}
real.T = summary(lm(Fitness ~ Delta.promoter.strength, data = promoter_yfp_data))$coefficients[6]
# real.T = summary(lm(Fitness ~ log(Delta.promoter.strength), data = promoter_yfp_data))$coefficients[6]
(P.value = sum(real.T<T.dist)/N) 
# P = 0.1121 (linear)
# P = 0.0882 (log)

plot(log(promoter_yfp_data$Delta.promoter.strength), (promoter_yfp_data$Fitness))
summary(mod<-lm(Fitness~log(Delta.promoter.strength), data = promoter_yfp_data))
abline(mod)


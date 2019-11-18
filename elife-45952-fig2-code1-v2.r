# Figure 2C: YFP vs fitness

fig2c = read_excel("Fig 2C - YFP-fitness data.xlsx", 
            +     skip = 1)

fig2c$Mutation = as.factor(fig2c$Mutation)

lm.out = lm(Mean.YFP ~ Mean.fitness, data = fig2c)


newx = seq(min(fig2c$Mean.fitness), max(fig2c$Mean.fitness),
           length.out = 100)
preds = predict(lm.out, newdata = data.frame(Mean.fitness = newx),
                interval = 'confidence')
library(scales)
myGrey = rgb(220,220,220, maxColorValue = 255)

plot(Mean.YFP ~ Mean.fitness, data = fig2c,
     ylab = "Relative Maximum Fluorescence", xlab = "Relative Fitness",
     pch = 18)
polygon(c(rev(newx), newx),
        c(rev(preds[,3]), preds[,2]),
        col= myGrey, border = NA)
points(fig2c$Mean.fitness, fig2c$Mean.YFP, pch = 18)
abline(v = 1, lty = 2)
abline(h = 1, lty = 2)
abline(lm.out, col = 'blue', lwd = 2)


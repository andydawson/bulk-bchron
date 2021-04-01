library(plyr)

cc = read.csv('data/chroncontrol_summary_pollen_full.csv')

ncontrols = ddply(cc, .(datasetid), nrow)

id_keep = ncontrols$datasetid[which(ncontrols$V1>1)]

cc = cc[which(cc$datasetid %in% id_keep),]

##
## core tops
##

core_tops = cc[which(cc$type %in% c('Core top', 'Collection date')),]
core_tops$error = abs(core_tops$limitolder - core_tops$limityounger)
foo = core_tops[,c('datasetid', 'siteid', 'sitename', 'age', 'depth', 'limityounger', 'limitolder', 'error'),]

foo = foo[which(foo$error > 100),]
write.csv(foo, 'data/core_tops.csv', row.names=FALSE)

ggplot() + geom_point(data=core_tops, aes(x=age, y=error))


hist(core_tops$error, breaks=30)

ct_error = mean(core_tops$error, na.rm=TRUE)

##
## tephra
##

tephra = cc[which(cc$type %in% c('Tephra')),]
tephra$error = abs(tephra$limitolder - tephra$limityounger)
foo = tephra[,c('datasetid', 'siteid', 'sitename', 'age', 'depth', 'limityounger', 'limitolder', 'error'),]

# foo = foo[which(foo$error > 100),]
# write.csv(foo, 'data/core_tops.csv', row.names=FALSE)

ggplot() + geom_point(data=tephra, aes(x=age, y=error))


hist(tephra$error, breaks=30)

tephra_error = mean(tephra$error, na.rm=TRUE)

##
## radiocarbon dates
##

radio = cc[which(cc$type %in% c('Radiocarbon')),]#, 'Radiocarbon, reservoir correction', 'Radiocarbon, average of two or more dates')),]
radio$error = abs(radio$limitolder - radio$limityounger)/2

radio = radio[which((!is.na(radio$error)) & (radio$error>0)),]

write.csv(radio, 'data/radiocarbon-dates-errors.csv', row.names=FALSE)


ggplot() + 
  geom_point(data=radio, aes(x=age, y=error)) + 
  geom_smooth(data=radio, aes(x=age, y=error)) +
  theme_bw() +
  xlab('Radiocarbon age (YBP)') + 
  ylab('Standard deviation') + 
  ylim(c(0,7500))

library(mgcv)

hist(radio$error, breaks=50, xlim=c(0,5000))

mod <- gam(error ~ s(age, k=15), data=radio, method='REML', family=Gamma(link="identity"))
plot(mod, shade = TRUE, seWithMean = TRUE, residuals = TRUE, pch = 16, cex = 0.8)
gam.check(mod)


library(gratia)
age_new = seq(-20, 52500, length.out = 1000)




draw(mod)

new_age = data.frame(age=age_new)

model_p = cbind(new_age, data.frame(predict.gam(mod, new_age, type='response', se.fit=TRUE)))
crit.t = qt(0.975, df = mod$df.residual)
new_age = transform(model_p,
                      upper = fit + (crit.t * se.fit),
                      lower = fit - (crit.t * se.fit))

ggplot() +
  geom_point(data=radio, aes(x=age, y=error)) + 
  geom_line(data=new_age, aes(x=age, y=fit), colour='blue') +
  geom_ribbon(data=new_age, aes(x=age, ymin=lower, ymax=upper), fill="grey", alpha=0.5) + 
  theme_bw() +
  xlab('Radiocarbon age (YBP)') + 
  ylab('Standard deviation')


nsim = 20
sims <- simulate(mod, nsim = nsim, newdata = new_age, unconditional = TRUE)
## rearrange the output into a long/tidy format
colnames(sims) <- paste0("sim", seq_len(nsim))
sims <- setNames(stack(as.data.frame(sims)), c("simulated", "run"))
sims <- transform(sims, age = rep(new_age$age, nsim),
                  simulated = simulated)
## Plot simulated trends
smallSim.plt <- ggplot(new_age, aes(x = age, y = fit)) +
  geom_line(data = sims,
            mapping = aes(y = simulated, x = age, group = run),
            colour = "grey80") +
  geom_line(lwd = 1) +
  labs(y = 'SD', x = "Radiocarbon age (YBP)") + 
  theme_bw()
smallSim.plt




# ##
# 
# radio = radio[order(radio$age),]
# 
# 
# mod <- gamm(error ~ s(age), 
#             data=radio,
#             # correlation = corAR1(form = ~age),
#             method='REML',
#             family=Gamma(link="identity"),
#             control = list(niterEM = 0, optimMethod = "BFGS", opt="optim"))
# 
# plot(mod, shade = TRUE, seWithMean = TRUE, residuals = TRUE, pch = 16, cex = 0.8)
# gam.check(mod$gam)
# 
# 
# library(gratia)
# age_new = seq(-20, 52500, length.out = 1000)
# 
# new_age = data.frame(age=age_new)
# 
# model_p = cbind(new_age, data.frame(predict(mod$gam, new_age, se.fit=TRUE)))
# crit.t = qt(0.975, df = df.residual(mod$gam))
# new_age = transform(model_p,
#                     upper = fit + (crit.t * se.fit),
#                     lower = fit - (crit.t * se.fit))
# 
# foo = acf(model_p$fit, lag.max=500)
# 
# 
# ggplot() +
#   geom_point(data=radio, aes(x=age, y=error)) +
#   geom_line(data=new_age, aes(x=age, y=fit), colour='blue') +
#   geom_ribbon(data=new_age, aes(x=age, ymin=lower, ymax=upper), fill="grey", alpha=0.5) +
#   theme_bw() +
#   xlab('Radiocarbon age (YBP)') +
#   ylab('Standard deviation')

# 
# nsim = 20
# sims <- simulate(mod, nsim = nsim, newdata = new_age, unconditional = TRUE)
# ## rearrange the output into a long/tidy format
# colnames(sims) <- paste0("sim", seq_len(nsim))
# sims <- setNames(stack(as.data.frame(sims)), c("simulated", "run"))
# sims <- transform(sims, age = rep(new_age$age, nsim),
#                   simulated = simulated)
# ## Plot simulated trends
# smallSim.plt <- ggplot(new_age, aes(x = age, y = fit)) +
#   geom_line(data = sims,
#             mapping = aes(y = simulated, x = age, group = run),
#             colour = "grey80") +
#   geom_line(lwd = 1) +
#   labs(y = 'SD', x = "Radiocarbon age (YBP)") + 
#   theme_bw()
# smallSim.plt
# preds = data.frame(age=new_age, preds=model_p$fit, se=model_p$se.fit)
# ggplot() +
#   geom_point(data=radio, aes(x=age, y=error)) + 
#   geom_line(data=preds, aes(x=age, y=preds), colour='blue') +
#   geom_ribbon(data=preds, aes(x=age, ymin=preds - 2*se, ymax=preds + 2*se), fill="grey", alpha=0.5) + 
#   theme_bw() +
#   xlab('Radiocarbon age (YBP)') + 
#   ylab('Standard deviation')
# 
# 
# ggplot() +
#   geom_point(data=radio, aes(x=age, y=error)) + 
#   geom_line(data=new_age, aes(x=age, y=fit), colour='blue') +
#   geom_ribbon(data=new_age, aes(x=age, ymin=lower, ymax=upper), fill="grey", alpha=0.5) + 
#   theme_bw() +
#   xlab('Radiocarbon age (YBP)') + 
#   ylab('Standard deviation')

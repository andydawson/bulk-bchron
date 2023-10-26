library(plyr)
library(ggplot2)
library(mgcv)
library(gratia)
library(reshape2)
library(dplyr)

cc = read.csv('data/chroncontrol_summary_pollen_full.csv')

# number of controls for each dataset
ncontrols = ddply(cc, .(datasetid), nrow)

# for error analysis 
# only keep datasets with more than one control 
id_keep = ncontrols$datasetid[which(ncontrols$V1>1)]

cc = cc[which(cc$datasetid %in% id_keep),]

# fix data entry errors
cc[which(cc$chroncontrolid == 105848),'age'] = -10
cc[which(cc$chroncontrolid == 105849),'age'] = 12
cc[which(cc$chroncontrolid == 105848),'limityounger'] = 1950 - cc[which(cc$chroncontrolid == 105848),'limityounger']
cc[which(cc$chroncontrolid == 105849),'limityounger'] = 1950 - cc[which(cc$chroncontrolid == 105849),'limityounger']

cc[which(cc$chroncontrolid == 104870),'type'] = 'Tephra'

cc[which(cc$limityounger < (-70)), 'limityounger'] = -70
cc[which((cc$type=="Lead-210")&(cc$limitolder>200)), 'limitolder'] = 200

cc[which(is.na(cc$age)),'age'] = (cc[which(is.na(cc$age)),'limityounger'] + cc[which(is.na(cc$age)),'limitolder']) / 2

cc = cc[which(!is.na(cc$age)),]

######################################################################################################
## summarize controls
######################################################################################################



######################################################################################################
## core top error
######################################################################################################

core_tops = cc[which(cc$type %in% c('Core top', 'Collection date')),]
core_tops$error = abs(core_tops$limitolder - core_tops$limityounger)
core_tops = core_tops[,c('datasetid', 'siteid', 'sitename', 'age', 'depth', 'limityounger', 'limitolder', 'error'),]

core_tops_problem = core_tops[which(core_tops$error > 100),]

ggplot() + geom_point(data=core_tops, aes(x=age, y=error))


hist(core_tops$error, breaks=30)

ct_error = mean(core_tops$error, na.rm=TRUE)
ct_error

######################################################################################################
## tephra
######################################################################################################

tephra = cc[which(cc$type %in% c('Tephra')),]
tephra$error = abs(tephra$limitolder - tephra$limityounger)
foo = tephra[,c('datasetid', 'siteid', 'sitename', 'age', 'depth', 'limityounger', 'limitolder', 'error'),]

# foo = foo[which(foo$error > 100),]
# write.csv(foo, 'data/core_tops.csv', row.names=FALSE)

ggplot() + geom_point(data=tephra, aes(x=age, y=error))
ggplot() + geom_point(data=subset(tephra, error<5000), aes(x=age, y=error))

hist(tephra$error, breaks=30)

tephra_error = mean(tephra$error, na.rm=TRUE)
tephra_error


######################################################################################################
## lead 210
######################################################################################################

lead = cc[which(cc$type %in% c('Lead-210')),]
lead$error = abs(lead$limitolder - lead$limityounger)

lead = lead[which((!is.na(lead$error)) & (lead$error>0)),]

lead$age_positive = 1950 - lead$age 

write.csv(lead, 'data/lead-dates-errors.csv', row.names=FALSE)

ggplot() + 
  geom_point(data=lead, aes(x=age, y=error)) + 
  geom_smooth(data=lead, aes(x=age, y=error)) +
  theme_bw() +
  xlab('Lead-210 Age') + 
  ylab('Standard deviation') + 
  ylim(c(0,300))

mod <- gam(error ~ s(age, k=4), data=lead, method='REML', family=Gamma(link="log"))
# mod <- gam(error ~ s(age_positive, k=4), data=lead, method='REML', family=Gamma(link="log"))
plot(mod, shade = TRUE, seWithMean = TRUE, residuals = TRUE, pch = 16, cex = 0.8)
gam.check(mod)

age_new = seq(-66, 206, length.out = 1000)

new_age = data.frame(age=age_new)

# model_p = cbind(new_age, data.frame(predict.gam(mod, new_age, type='response', se.fit=TRUE)))
# crit.t = qt(0.975, df = mod$df.residual)
# new_age = transform(model_p,
#                     upper = fit + (crit.t * se.fit),
#                     lower = 0)#ifelse((fit - (crit.t * se.fit))>0, fit - (crit.t * se.fit), 0))

ggplot() +
  geom_point(data=lead, aes(x=age, y=error)) + 
  geom_line(data=new_age, aes(x=age, y=fit), colour='blue') +
  geom_ribbon(data=new_age, aes(x=age, ymin=lower, ymax=upper), fill="grey", alpha=0.5) + 
  theme_bw() +
  xlab('Radiocarbon age (YBP)') + 
  ylab('Standard deviation')

lead_sims = simulate(mod, nsim = 1, data = new_age)

lead_sims_df = data.frame(age = new_age[,'age'], lead_sims)
foo_melt = melt(lead_sims_df, id.vars = c('mar'))

cal_sim_gam_sum = foo_melt %>% 
  group_by(mar) %>%
  summarize(alb_mean = mean(value), 
            alb_lo = quantile(value, c(0.025)), 
            alb_mid = quantile(value, c(0.5)), 
            alb_hi = quantile(value, c(0.975)))

ggplot() +
  geom_point(data=lead, aes(x=age, y=error)) + 
  geom_line(data=new_age, aes(x=age, y=fit), colour='blue') +
  geom_ribbon(data=new_age, aes(x=age, ymin=lower, ymax=upper), fill="grey", alpha=0.5) + 
  theme_bw() +
  xlab('Radiocarbon age (YBP)') + 
  ylab('Standard deviation')

######################################################################################################
## radiocarbon
######################################################################################################

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

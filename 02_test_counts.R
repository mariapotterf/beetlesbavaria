
# test data


load(file = "outData/ips.Rdata")



# ANOVA : # does the DOY of max increase differes between years?  ----------------
#http://www.sthda.com/english/wiki/one-way-anova-test-in-r#check-your-data
# 
# one way anova - finds that there are differenes, we don't know where
res.aov <- aov(doy ~ factor(year), max.diff.doy.sf)

summary(res.aov)

# multiple way anova = pairwise comparison between groups
TukeyHSD(res.aov) # challengin to read differences, if many groups

pairwise.t.test(max.diff.doy.sf$doy, factor(max.diff.doy.sf$year),
                p.adjust.method = "BH")


### check ANOVA assumptions:
# 1. Homogeneity of variances
# 2. Normality


# 1. Homogeneity of variances
plot(res.aov, 1)

# test for homonenity f variance: 
library(car)
# high p-value = variances are the same (H0) 
leveneTest(doy ~ factor(year), data = max.diff.doy.sf)
# low p-values, so variancs are different ->

# use ANOVA with no assumption of equal variances
res.aov_no_equal <- oneway.test(doy ~ factor(year), data = max.diff.doy.sf)



### t-test for not equal variances : ---------------------------
pairwise.t.test(max.diff.doy.sf$doy, factor(max.diff.doy.sf$year),
                p.adjust.method = "BH", pool.sd = FALSE)


# 2. Normality
plot(res.aov, 2)

# check by Shapiro-Wilk test for residuals
# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test # does not work, have more values then 5000
shapiro.test(x = aov_residuals )  # not normaly distributed


# run rather Kruskal-Wallis test - not equal variances
kruskal.test(doy ~ factor(year), data = max.diff.doy.sf)


# post-hoc test: pairwise-wilcox
wilk <-pairwise.wilcox.test(max.diff.doy.sf$doy, factor(max.diff.doy.sf$year),
                            p.adjust.method = "BH")
# get matrix of p-values for further ploting
wilk$p.value

p.vals <- replace(wilk$p.value,wilk$p.value>0.05,1)


# show p-values - first transformation is needed, 
# as this is the opposite of the corr() plot - larger is better,
# while in p-values smaller is better
# https://stackoverflow.com/questions/61861049/how-plot-results-from-a-pairwise-wilcox-test
my_transformed_pvals=-log10(p.vals)  # negative of 'log' works, as all p-vals are <1

library ( corrplot)
corrplot(as.matrix(my_transformed_pvals),is.corr=F)
corrplot(as.matrix(p.vals),is.corr=F)


# !!! my groups are not independent! they are colected on teh same traps, 
# they are autocorrelated!
# use anova with repeated measures:

aov_rep <- aov(doy~factor(year)+Error(factor(globalid)), data = max.diff.doy.sf)
summary(aov_rep)



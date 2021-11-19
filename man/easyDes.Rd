\name{easyDes}
\alias{easyDes}
\title{An Easy Way to Descriptive Analysis}
\description{
  Descriptive analysis is essential for publishing medical articles.
  This package provides an easy way to conduct the descriptive analysis.
  1. Both numeric and factor variables can be handled. For numeric variables, normality test will be applied to choose the parametric and nonparametric test.
  2. Both two or more groups can be handled. For groups more than two, the post hoc test will be applied, 'Tukey' for the numeric variables and 'FDR' for the factor variables.
  3. T test, ANOVA or Fisher test can be forced to apply.
  4. Mean and standard deviation can be forced to display.
}
\keyword{Descriptive analysis}
\usage{
easyDes(nc.g,nc.n,nc.f,nc.of,dataIn,fisher,aov,t,mean,
    mcp.test.method,mcp.stat,mcp.t.test,mcp.t.test.method,
    table.margin,decimal.p,decimal.prop)
}
\arguments{
  \item{nc.g}{integer, the column number of the grouping variable, length of 'nc.g' must be 1}
  \item{nc.n}{numeric vector, the column number of the numeric variable, length of 'nc.n' can be more than 1}
  \item{nc.f}{numeric vector, the column number of the factor variable, length of 'nc.f' can be more than 1}
  \item{nc.of}{numeric vector, the column number of the ordinal factor variable, length of 'nc.of' can be more than 1}
  \item{dataIn}{data frame including variables above}
  \item{fisher}{logic, whether to apply Fisher test by force, the default is 'TRUE'}
  \item{aov}{logic, whether to apply ANOVA test by force, the default is 'FALSE'}
  \item{t}{logic, whether to apply T test by force, the default is 'FALSE'}
  \item{mean}{logic, whether to disply the mean and standar deviation for the numeric variables by force, the default is 'FALSE'}
  \item{mcp.test.method}{character, specific for ANOVA, the method for the multiple comparisons in 'multcomp' package, 'Tukey' or 'Dunnett'}
  \item{mcp.stat}{logic, whether to display the statistic for the multiple comparsions}
  \item{mcp.t.test}{logic, specific for ANOVA, wether to use the pairwise t tests as the multiple comparsions instead of that in 'multcomp' package}
  \item{mcp.t.test.method}{character, specific for 'mcp.t.test'==TRUE, the method for the pairwise t tests, 'holm' (Holm, 1979), 'hochberg' (Hochberg, 1988), 'hommel' (Hommel, 1988), 'bonferroni', 'BH' (Benjamini & Hochberg, 1995), 'BY' (Benjamini & Yekutieli, 2001), 'fdr', 'none'}
  \item{table.margin}{1 or 2, which margin of the table should be calculated the proportion, 1=row, 2=column}
  \item{decimal.p}{integer, the number of decimals of the p value}
  \item{decimal.prop}{integer, the number of decimals of the proportions for factor variables}
}
\value{
  \item{total}{the descriptive statistic for all data}
  \item{group names}{the descriptive statistic for data of each group}
  \item{method}{the method applied to test between groups, i.e. ANOVA and Tukey, Fisher and FDR}
  \item{statistic}{the statistic of test, i.e., the 'W' to Wilcoxon test, the 'chi-squared' to Kruskal-Wallis, the 't' to t test, the 'F' to ANOVA test}
  \item{p.value}{the p value derived from the test between groups}
  \item{stat.*_va_*}{the statistic derived from the post hoc test, the 't' of Tukey for ANOVA, the 'q' of Nemenyi for Kruskal-Wallis}
  \item{p.*_va_*}{the p value derived from the post hoc test}

}
\author{
Zhicheng Du<dgdzc@hotmail.com>, Yuantao Hao<haoyt@mail.sysu.edu.cn>
}
\details{
1. Nemenyi test was used as a Kruskal-Wallis post-hoc test.

2. FDR (False Discovery Rate) was used to adjust the p values after pairwise comparision of Chi-square test or Fisher test.

3. Tukey test was used as a ANOVA (Analysis of Variance) post-hoc test.

4. Shapiro-Wilk test was used as normality test if the sample size was between 3~5,000, while Kolmogorov-Smirnov test was used if the sample size was greater than 5,000.
}
\note{
Please feel free to contact us, if you have any advice and find any bug!

Update description:

Version 2.0: 1. T test can be forced to apply.

Version 3.0: 1. Fixing the wrong colnames in Chi-squre test. 2. Limiting the number of the decimal digits of the statistic in Chi-squre test to three. 3. The number of decimal digits of the propotion for the factor variables can be set free.

Version 4.0: 1. Mean and standard deviation can be forced to display. 2. The help document has been revised. 3. Fix the problem with more than 5,000 samples in the normality test.

Version 5.0: 1. Unify the number of decimal digits (i.e., output "0.010" rather than "0.01" for p value). 2. Add the 'nc.of' to analyze ordinal factors.

Version 6.0: 1. Add the pairwise t tests for the the multiple comparsions. 2. Fix the error of "Increase workspace or consider using 'simulate.p.value=TRUE'" in fisher test. 3. Add the 'table.margin' argument. 4. Add the 'decimal.p' argment. 5. Fix the bugs caused by the names with specific characters in numeric variables.

}
\examples{
group=rep(c(0,1),each=30)
nx1=rnorm(60)
nx2=rnorm(60)
fx1=rep(c(1:3),20)
fx2=rep(c(1:5),12)
fx3=factor(fx2)
data=data.frame(group,nx1,nx2,fx1,fx2,fx3)

easyDes(nc.g=1,nc.n=2:3,nc.f=4:5,nc.of=6,dataIn=data,
        fisher=TRUE,aov=FALSE,t=FALSE,mean=FALSE,mcp.stat=FALSE)
easyDes(nc.g=4,nc.n=2:3,nc.f=c(5,5),nc.of=6,dataIn=data,
        fisher=TRUE,aov=FALSE,t=FALSE,mean=FALSE,mcp.stat=FALSE)
easyDes(nc.g=4,nc.n=3,nc.f=5,nc.of=6,dataIn=data,
        fisher=TRUE,aov=FALSE,t=FALSE,mean=FALSE,mcp.stat=TRUE)
}


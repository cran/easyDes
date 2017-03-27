\name{easyDes}
\alias{easyDes}
\title{An Easy Way to Descriptive Analysis}
\description{
  Descriptive analysis is essential for publishing medical articles.
  This package provides an easy way to conduct the descriptive analysis.
  1. Both numeric and factor variables can be handled. For numeric variables, normality test will be applied to choose the parametric and nonparametric test.
  2. Both two or more groups can be handled. For groups more than two, the post hoc test will be applied, 'Tukey' for the numeric variables and 'FDR' for the factor variables.
  3. T test, ANOVA or Fisher test can be forced to apply.
}
\keyword{Descriptive analysis}
\usage{
easyDes(nc.g,nc.n,nc.f,dataIn,fisher,aov,t)
}
\arguments{
  \item{nc.g}{integer, the column number of the grouping variable, length of 'nc.g' must be 1}
  \item{nc.n}{numeric vector, the column number of the numeric variable, length of 'nc.n' can be more than 1}
  \item{nc.f}{numeric vector, the column number of the factor variable, length of 'nc.n' can be more than 1}
  \item{dataIn}{data frame including variables above}
  \item{fisher}{logic, whether to apply Fisher test by force}
  \item{aov}{logic, whether to apply ANOVA test by force}
  \item{t}{logic, whether to apply T test by force}
}
\value{
  \item{total}{the descriptive statistic for all data}
  \item{group names}{the descriptive statistic for data of each group}
  \item{method}{the method applied to test between groups}
  \item{p.value}{the p value derived from the test between groups}
  \item{*_va_*}{the p value derived from the post hoc test}
}
\author{
Zhicheng Du<dgdzc@hotmail.com>, Yuantao Hao<haoyt@mail.sysu.edu.cn>
}
\note{
Please feel free to contact us, if you have any advice and find any bug!

Update description:
Version 2.0: 1. T test can be forced to apply.
}
\seealso{
nothing
}
\examples{
group=rep(c(0,1),each=30)
nx1=rnorm(60)
nx2=rnorm(60)
fx1=rep(c(1:3),20)
fx2=rep(c(1:5),12)
data=data.frame(group,nx1,nx2,fx1,fx2)

easyDes(nc.g=1,nc.n=2:3,nc.f=4:5,dataIn=data,fisher=TRUE,aov=FALSE)
easyDes(nc.g=4,nc.n=2:3,nc.f=c(5,5),dataIn=data,fisher=TRUE,aov=FALSE)
}


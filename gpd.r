list.of.packages <- c("goft")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(goft)


argv = commandArgs(trailingOnly = TRUE)

y = scan(argv[1])
mexc = length(y)
p = gp_test(y)$p.value

mexc.c = seq(0,mexc-10,10)
z = y
i = 0

re = c()
for(i in mexc.c) {
  z = y[1:(mexc-i)]
  p = gp_test(z)$p.value
  re = rbind(re, c(i, p))
  if(!is.na(p) & p > 0.05) break
  i = i + 10
}

if(nrow(re) >=2) {
  mexc.c2 = seq(re[(nrow(re)-1),1]+1, re[nrow(re),1],1)
  re = c()
  for(i in mexc.c2) {
    z = y[1:(mexc-i)]
    p = gp_test(z)$p.value
    re = rbind(re, c(i, p))
    if(!is.na(p) & p > 0.05) break
    i = i + 1
  }
}

p = re[nrow(re),2]
len = mexc-re[nrow(re),1]

cat(c(p, len))



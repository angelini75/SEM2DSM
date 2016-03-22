# > summary(tk)$sigma
# [1] 0.9808548
# > summary(tk)$sigma
# [1] 0.9808548
# > summary(oc)$sigma
# [1] 0.8790718
# > summary(tb)$sigma
# [1] 0.9542115
# > summary(sa)$sigma
# [1] 0.9720069
# > summary(ea)$sigma
# [1] 0.9013179
# > summary(eb)$sigma
# [1] 0.8719843
# > summary(bt)$sigma
# [1] 0.8582147

setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/1_Sampling/Data")

mean((read.csv("res.thick.csv")[,10]^2)/(0.9808548*7.47)^2)
median((read.csv("res.thick.csv")[,10]^2)/(0.9808548*7.47)^2)


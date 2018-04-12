# INPUT PARAMETERS
# _X_SOURCE, _LO_LIMIT_1
# raw data in csv using db structure, i.e., a single row contains:
# code of representative: voter_id, code of division: vote_event_id, vote option according to Popolo standar: option (i.e., one of "yes", "no", "abstain", "not voting", "absent"))
# first row is the header
# for example:
# "voter_id","vote_event_id","option"
# “Joe Europe”,”Division-007”,”yes”
X_source = _X_SOURCE
# lower limit to eliminate from calculations, e.g., .1; number
lo_limit = _LO_LIMIT_1

# reorder data; divisions x persons
# we may need to install and/or load some additional libraries
# install.packages("reshape2")
library("reshape2")

X_source$vote_event_id = as.factor(X_source$vote_event_id)
X_source$voter_id = as.factor(X_source$voter_id)

X_source$option_numeric = rep(0,length(X_source$option))
X_source$option_numeric[X_source$option=='yes'] = 1
X_source$option_numeric[X_source$option=='no'] = -1
X_source$option_numeric[X_source$option=='abstain'] = -1    #may be 0 in some parliaments
X_source$option_numeric[X_source$option=='not voting'] = NA
X_source$option_numeric[X_source$option=='absent'] = NA
X_source$option_numeric = as.numeric(X_source$option_numeric)

#prevent reordering, which is behaviour of acast:
X_source$voter_id = factor(X_source$voter_id, levels=unique(X_source$voter_id))
X_raw = t(acast(X_source,voter_id~vote_event_id,fun.aggregate=mean,value.var='option_numeric'))
X_people = dimnames(X_raw)[[2]]
X_vote_events = dimnames(X_raw)[[1]]
#X_raw = apply(X_raw,1,as.numeric)
mode(X_raw) = 'numeric'



# WEIGHTS
# weights 1 for divisions, based on number of persons in division
w1 = apply(abs(X_raw)==1,1,sum,na.rm=TRUE)/max(apply(abs(X_raw)==1,1,sum,na.rm=TRUE))
w1[is.na(w1)] = 0
# weights 2 for divisions, "100:100" vs. "195:5"
w2 = 1 - abs(apply(X_raw==1,1,sum,na.rm=TRUE) - apply(X_raw==-1,1,sum,na.rm=TRUE))/apply(!is.na(X_raw),1,sum)
w2[is.na(w2)] = 0
# total weights:
w = w1 * w2

# analytical charts for weights:
plot(w1)
plot(w2)
plot(w)


# MISSING DATA
# index of missing data; divisions x persons
I = X_raw
I[!is.na(X_raw)] = 1
I[is.na(X_raw)] = 0


# EXCLUSION OF REPRESENTATIVES WITH TOO FEW VOTES (WEIGHTED)
# weights for non missing data; division x persons
I_w = I*w
# sum of weights of divisions for each persons; vector of length “persons”
s = apply(I_w,2,sum)
person_w = s/sum(w)
# index of persons kept in calculation; vector of length “persons”
person_I = person_w > lo_limit

# cutted (omitted) persons with too few weighted votes; division x persons
X_c = X_raw[,person_I]
# scale data; divisions x persons (mean=0 and sd=1 for each division); scaled cutted persons with too few weighted votes; division x persons
X_c_scaled = t(scale(t(X_c),scale=TRUE))
# scaled with NA->0 and cutted persons with too few weighted votes; division x persons
X_c_scaled_0 = X_c_scaled
X_c_scaled_0[is.na(X_c_scaled_0)] = 0
# weighted scaled with NA->0 and cutted persons with too few weighted votes; division x persons
X = X_c_scaled_0 * sqrt(w)  # X is shortcut for X_c_scaled_0_w

# “X’X” MATRIX
# weighted X’X matrix with missing values substituted and excluded persons; persons x persons
C = t(X) %*% X

# DECOMPOSITION
# eigendecomposition
Xe=eigen(C)
# W (rotation values of persons)
V = Xe$vectors
# projected divisions into dimensions
Xy = X %*% V

# analytical charts of projection of divisions and lambdas
plot(Xy[,1],Xy[,2])
plot(sqrt(Xe$values[1:min(10,dim(Xy))]))

# lambda matrix
sigma = sqrt(Xe$values)
sigma[is.na(sigma)] = 0
lambda = diag(sigma)

# projection of persons into dimensions
X_proj = V %*% lambda
# unit-standardized projection of persons into dimensions
X_proj_unit = X_proj / sqrt(apply(X_proj^2,1,sum))

# analytical charts
plot(X_proj[,1],X_proj[,2])
plot(X_proj_unit[,1],X_proj_unit[,2])

# lambda^-1 matrix
lambda_1 = diag(sqrt(1/Xe$values))
lambda_1[is.na(lambda_1)] = 0

# U (rotation values of divisions)
U = X %*% V %*% lambda_1

# analytical charts
# second projection
X_proj2 = t(X) %*% U
# second unit scaled projection of persons into dimensions
X_proj2_unit = X_proj2 / sqrt(apply(X_proj2^2,1,sum))
# they should be equal:
plot(X_proj[,1],X_proj2[,1])
#plot(X_proj[,2],X_proj2[,2])

# save first two dimensions with persons' ids:
#write.csv(cbind(X_people[person_I],X_proj_unit[,1:2]),file="output.csv")


# PROJECTION OF PART OF THE LEGISLATIVE TERM INTO THE WHOLE TERM MODEL
# additional parameters:
# _TI,   _LO_LIMIT_2
# lower limit to eliminate from projections (may be the same as lo_limit); number
lo_limit_T = _LO_LIMIT_2
# indices whether divisions are in the given part or not
# vector of length div, contains TRUE or FALSE
TI = _TI
# missing values treatment (from overall model)
X_raw_T_c = X_raw[,person_I]
X_raw_T_c[!TI,] = NA

# Indices of NAs; division x person
TI_c = X_raw_T_c
TI_c[!is.na(TI_c)] = 1
TI_c[is.na(TI_c)] = 0

# weights for non missing data; divisions x persons
TI_c_w = TI_c * w
# sum of weights of divisions for each persons
s = apply(TI_c_w,2,sum)
person_T_w = s/max(s)
# index of persons in calculation
person_TI = person_T_w > lo_limit_T

# scaling
XT_c_scaled = (X_raw_T_c-attr(X_c_scaled_0,"scaled:center"))/attr(X_c_scaled_0,"scaled:scale")
# weighting
XT_c_scaled_w_0 = XT_c_scaled * sqrt(w)
XT_c_scaled_w_0[is.na(XT_c_scaled_w_0)] = 0

# weighted scaled with NA->0 and cutted persons with too few votes division x persons
XT = XT_c_scaled_w_0[,person_TI]    # ~XT_cc_scaled_w_0
# index of missing data cutted persons with too few votes divisions x persons
TI_cc = TI_c[,person_TI]

# weights of divisions for each person, person x division
X_cc = X[,person_TI]
aU = abs(U)
aX_cc = abs(X_cc)
aX_proj_cc = t(aX_cc) %*% aU
aXT_cc_scaled_w_0 = abs(XT)
aXT_proj_cc = t(aXT_cc_scaled_w_0) %*% aU
div_person_weights = aXT_proj_cc/aX_proj_cc

# projection of persons
XT_proj = t(XT) %*% U / div_person_weights
# standardization (unit)
XT_proj_unit = XT_proj / sqrt(apply(XT_proj^2,1,sum))

# analytical chart:
plot(XT_proj_unit[,1],XT_proj_unit[,2])


# PROJECTION OF ADDITIONAL PERSONS INTO THE WHOLE TERM MODEL
# additional persons:
# _XAP_SOURCE
# the same structure as X_source
# needs to have the same vote events (and their order!) as the main analysis!
Xap_source = _XAP_SOURCE

# reorder, as with X_source
Xap_source$vote_event_id = as.factor(Xap_source$vote_event_id)
Xap_source$voter_id = as.factor(Xap_source$voter_id)

Xap_source$option_numeric = rep(0,length(Xap_source$option))
Xap_source$option_numeric[Xap_source$option=='yes'] = 1
Xap_source$option_numeric[Xap_source$option=='no'] = -1
Xap_source$option_numeric[Xap_source$option=='abstain'] = -1    #may be 0 in some parliaments
Xap_source$option_numeric[Xap_source$option=='not voting'] = NA
Xap_source$option_numeric[Xap_source$option=='absent'] = NA
Xap_source$option_numeric = as.numeric(Xap_source$option_numeric)

#prevent reordering, which is behaviour of acast:
Xap_source$voter_id = factor(Xap_source$voter_id, levels=unique(Xap_source$voter_id))
Xap_raw = t(acast(Xap_source,voter_id~vote_event_id,fun.aggregate=mean,value.var='option_numeric'))
Xap_people = dimnames(Xap_raw)[[2]]
Xap_vote_events = dimnames(Xap_raw)[[1]]
#Xap_raw = apply(Xap_raw,1,as.numeric)
mode(Xap_raw) = 'numeric'

# standardize using parameters from previous standardization:
Xap_scaled = (Xap_raw-attr(X_c_scaled_0,"scaled:center"))/attr(X_c_scaled_0,"scaled:scale")
Xap_scaled[is.na(Xap_scaled)] = 0
Xap_scaled[is.nan(Xap_scaled)] = 0
Xap_scaled[Xap_scaled==-Inf] = 0
Xap_scaled[Xap_scaled==Inf] = 0

Xap_scaled_w = Xap_scaled * sqrt(w)

# projection of additional persons into dimensions, 2nd type:
Xap_proj2 = t(Xap_scaled_w) %*% U
# standardized projection of additional persons into dimensions, 2nd type:
Xap_proj2_unit = Xap_proj2 / sqrt(apply(Xap_proj2^2,1,sum))

# analytical charts
plot(Xap_proj2[,1],Xap_proj2[,2])
plot(Xap_proj2_unit[,1],Xap_proj2_unit[,2])


# CUTTING LINES
# additional parameters:
# _N_FIRST_DIMENSIONS
# how many dimensions?
n_first_dimensions = _N_FIRST_DIMENSIONS

# loss function
LF = function(beta0) -1*sum(apply(cbind(y*(x%*%beta+beta0),zeros),1,min))

# preparing variables
normals = Xy[,1:n_first_dimensions]
loss_f = data.frame(matrix(0,nrow=dim(X_raw)[1],ncol=4))
colnames(loss_f)=c("Parameter1","Loss1","Parameter_1","Loss_1")
parameters = data.frame(matrix(0,nrow=dim(X_raw)[1],ncol=3))
colnames(parameters)=c("Parameter","Loss","Direction")

# x-values
#xfull = t(t(Xe$vectors[,1:n_first_dimensions]) * sqrt(Xe$values[1:n_first_dimensions]))
#xfull = X_proj[,1:n_first_dimensions]
# unit x-values
xfullu = X_proj_unit[,1:n_first_dimensions]


#calculating all cutting lines
for (i in as.numeric(1:dim(X_raw)[1])) {
  beta = Xy[i,1:n_first_dimensions]
  y = t(as.matrix(X_raw[i,]))[,person_I]
  x = xfullu[which(!is.na(y)),]
  y = y[which(!is.na(y))]
  zeros = as.matrix(rep(0,length(y)))
  # note: “10000” should be enough for any real-life case:
  res1 = optim(c(1),LF,method="Brent",lower=-10000,upper=10000)        
  # note: the sign is arbitrary, the real result may be -1*; we need to minimize the other way as well
  y=-y
  res2 = optim(c(1),LF,method="Brent",lower=-10000,upper=10000) 
  
  # the real parameter is the one with lower loss function
  # note: theoretically should be the same (either +1 or -1) for all divisions(rows), however, due to missing data, the calculation may lead to a few divisions with the other direction
  loss_f[i,] = c(res1$par,res1$value,res2$par,res2$value)
  if (res1$value<=res2$value) {
    parameters[i,] = c(res1$par,res1$value,1)
  } else {
    parameters[i,] = c(res2$par,res2$value,-1)
  }
}
CuttingLines = list(normals=normals,parameters=parameters,loss_function=loss_f,weights=cbind(w1,w2))

# analytical charts
# cutting lines
# additional parameters:
# _MODULO, _LO_LIMIT_3
# to show only each _MODULO-division (for huge numbers of divisions may be useful to use it) # if set to 1, every division is shown
# _LO_LIMIT_3 is a lower limit used to plot only important divisions; number between [0,1]
lo_limit3 = _LO_LIMIT_3
modulo = _MODULO
plot(X_proj_unit[,1],X_proj_unit[,2])
I = w1*w2 > lo_limit3
for (i in as.numeric(1:dim(X_raw)[1])) {
  if (I[i] && ((i %% modulo) == 0)) {
    beta = CuttingLines$normals[i,]
    beta0 = CuttingLines$parameters$Parameter[i]
    abline(-beta0/beta[2],-beta[1]/beta[2])
  }
}
# normals of cutting lines, possibly with some limitations, e.g. 50, 20: 
plot(CuttingLines$parameters$Parameter/CuttingLines$normals[,1],ylim=c(-50,50))
plot(CuttingLines$parameters$Parameter/CuttingLines$normals[,2],ylim=c(-20,20))
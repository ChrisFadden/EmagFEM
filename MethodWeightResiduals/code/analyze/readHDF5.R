library(rhdf5)
DATAFP = "../../data/";

FP = paste0(DATAFP,"vector.h5");

#read the root group
mydata <- h5read(FP,"/")
print(mydata)

H5close()

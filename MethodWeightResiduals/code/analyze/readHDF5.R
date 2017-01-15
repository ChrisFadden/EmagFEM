library(rhdf5)
DATAFP = "../../data/";

FP = paste0(DATAFP,"vector.h5");

mydata <- h5read(FP,"/")

print(mydata$vecB)

H5close()

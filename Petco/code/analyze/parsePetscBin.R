#PARSE VECTOR
parsePetscBinVec <- function(filename,prefix="../../data/"){
  FP = paste0("../../data/",filename);
  to.read = file(FP, "rb")

  VecReadCorrect = 1211214
  #confirm vector read
  ierr  <- readBin(to.read,integer(),endian="big")
  stopifnot(ierr == VecReadCorrect)

  N <- readBin(to.read,integer(),endian="big")
  v <- readBin(to.read,double(),n=N,endian="big")
  return(v)
}

#PARSE MATRIX
parsePetscBinMat <- function(filename,prefix="../../data/"){
  FP = paste0(prefix,filename) 
  to.read = file(FP, "rb")

  MatReadCorrect = 1211216
  #confirm matrix read
  ierr  <- readBin(to.read,integer(),endian="big")
  stopifnot(ierr == MatReadCorrect)

  M <- readBin(to.read,integer(),endian="big")
  N <- readBin(to.read,integer(),endian="big")
  Total <- readBin(to.read,integer(),endian="big")

  rownz <- readBin(to.read,integer(),n=M,endian="big")
  colidx <- readBin(to.read,integer(),n=Total,endian="big")
  val <- readBin(to.read,double(),n=Total,endian="big")

  mat <- t(matrix(val,nrow=M,ncol=N))
  return(mat)
}


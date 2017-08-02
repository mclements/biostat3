require(epiR)
ir.est <- function(fail.str,py.str,group.str,data,method="exact"){
	sums <- by(data[c(fail.str,py.str)],data[[group.str]],colSums)
	dy <- matrix(unlist(sums),ncol=2,byrow=TRUE)
	colnames(dy) <- c("D","Y")
	rownames(dy) <- levels(factor(data[[group.str]]))
	cir <- epi.conf(dy, ctype="inc.rate", method)
	colnames(cir) <- c("IR","CI lower","CI upper")	
	cbind(dy,cir)
}
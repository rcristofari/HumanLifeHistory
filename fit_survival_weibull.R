
# Simulated dataset
a <- 5.55
b <- 87.56
x <- seq(1,100,1)
y <- jitter(weisv(x,a,b),5000)

# Subsample randomly:
# s <- sample((1:100), 15, replace=F)
# x <- x[s]
# y <- y[s]

fit_survival_weibull <- function(x,y, verbose=F, plot=F, converg=100, lsc=1e-6, usc=100, lsh=1e-6, ush=100, sbs=F){

# Weibull survival function:
weisv <- function(x,a,b){return(exp(-(x/b)^a))}

# Sum-of-squares from a Weibull survival function
weilsr <- function(x,y,a,b){
	sum((weisv(x,a,b)-y)^2)}

rowmin <- lsc
rowmax <- usc
colmin <- lsh
colmax <- ush
itermeanA <- vector()
itermeanB <- vector()
minima <- vector()
iter <- 1
criterion <- 0

while(criterion != 1){
		
# Build a 100 per 100 parameter matrix
lsr <- matrix(nrow=100, ncol=100, data=0)
as <- seq(rowmin,rowmax, length.out=100)
bs <- seq(colmin,colmax, length.out=100)

# Proceed step by step and pring the vectors along the way
if(sbs == T){
invisible(readline(prompt="Press [enter] to continue"))
print(as)
print(bs)
}

# Filling in the matrix with sum of squares
for(i in 1:100){
	for(j in 1:100){
		lsr[i,j]<-weilsr(x,y,as[i],bs[j])	
}}

# Identidying the LSS point and the associated parameter values
minimum <- which(lsr == min(lsr), arr.ind = TRUE)
minima[iter]<-min(lsr)

if(min(minimum[,1])>1){
rowmin <- as[min(minimum[,1])-1]
} else {
rowmin <- as[min(minimum[,1])]}
if(rowmin <= 0){rowmin<-0}

if(max(minimum[,1])<98){
	rowmax <- as[max(minimum[,1])+1]
	} else if (max(minimum[,1])>=98 && max(minimum[,1])<=100){
	rowmin <- 1e-6
	rowmax <- max(as)*2
	colmin <- 1e-6
	colmax <- max(bs)*2
	} else { rowmax <- as[max(minimum[,1])]}
if(rowmax <= 0){rowmax <-0}

if(min(minimum[,2])>1){
colmin <- bs[min(minimum[,2])-1]
} else {
colmin <- bs[min(minimum[,2])]}
if(colmin <= 0){colmin <-0}

if(max(minimum[,2])<98){
	colmax <- bs[max(minimum[,2])+1] 
	} else if (max(minimum[,2])>=98 && max(minimum[,2])<=100){
	colmin <- 1e-6
	colmax <- max(bs)*2
	rowmin <- 1e-6
	rowmax <- max(as)*2
	} else { colmax <- bs[max(minimum[,2])]}
if(colmax <= 0){colmax <-0}

if(sbs == T){
print(paste("IdA1", min(minimum[,1]), "IdA2",max(minimum[,1]),"IdB1" ,min(minimum[,2]),"IdB2", max(minimum[,2])), na.print="na")
print(paste("A-min:", rowmin, "A-max:",rowmax,"B-min:" ,colmin,"B-max:", colmax), na.print="na")}

itermeanA[iter] <- mean(c(as[min(minimum[,1])], as[max(minimum[,1])]))
itermeanB[iter] <- mean(c(bs[min(minimum[,2])], bs[min(minimum[,2])]))

# Convergence criterion (stability over 5 iterations and accuracy of min/max)
if(iter>6){
	if(( (itermeanA[iter] - mean(as) <= itermeanA[iter]/converg) 
		&& (itermeanB[iter] - mean(bs)) <= itermeanB[iter]/converg)){
				criterion <- 1 }}
		
# Details: cute plotting and verbosity level
#plot(raster(lsr))
if(plot==T){
	plot(y~x)
	curve(weisv(x, rowmin, colmin), from=min(x), to=max(x), col='red', add=T)}

if(verbose==T){
	print(paste("Iteration number ", iter))
	print(paste("Minimum least SS: ", min(lsr)))
	print(paste("Current parameter estimates: Shape ~", round(itermeanA[iter],2),' and Scale ~', round(itermeanB[iter],2), collapse='', sep=''))}

iter <- iter+1 } # End of the main While loop

if (plot ==T){
curve(weisv(x, a, b), from=min(x), to=max(x), col='orange', lty='dashed', add=T)}

return(data.frame(shape=itermeanA[iter-1], scale=itermeanB[iter-1]))


}


# Simulated dataset
a <- 2.3
b <- 37.56
x <- seq(1,100,1)
y <- jitter(weisv(x,a,b),1800)

# Subsample randomly:
s <- sample((1:100), 15, replace=F)
x <- x[s]
y <- y[s]

fit_survival_weibull <- function(x,y, verbose=F, plot=F, converg=100, lsc=1, usc=99, lsh=1, ush=99){

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
as <- seq(rowmin-(rowmin/10),rowmax+(rowmax/10), length.out=100)
bs <- seq(colmin-(colmin/10),colmax+(colmax/10), length.out=100)

# Filling in the matrix with sum of squares
for(i in 1:100){
	for(j in 1:100){
		lsr[i,j]<-weilsr(x,y,as[i],bs[j])	
}}

# Identidying the LSS point and the associated parameter values
minimum <- which(lsr == min(lsr), arr.ind = TRUE)
minima[iter]<-min(lsr)
rowmin <- as[min(minimum[,1])]
rowmax <- as[max(minimum[,1])]
colmin <- bs[min(minimum[,2])]
colmax <- bs[max(minimum[,2])]
itermeanA[iter] <- mean(c(rowmin, rowmax))
itermeanB[iter] <- mean(c(colmin, colmax))

# Convergence criterion (stability over 5 iterations and accuracy of min/max)
if(iter>6){
	if( ((itermeanA[iter] - mean(itermeanA[iter-5:iter])) <= itermeanA[iter]/converg) 
		&& ((itermeanB[iter] - mean(itermeanB[iter-5:iter])) <= itermeanB[iter]/converg) 
		&& ((rowmax-rowmin) < (rowmin/converg)) 
		&& ((colmax-colmin) < (colmin/converg))){
				criterion <- 1 }}
		
# Details: cute plotting and verbosity level
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
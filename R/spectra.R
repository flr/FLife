setGeneric('spectra', function(object,...)  standardGeneric('spectra'))

#Spectral analysis function
spectraFn<-function(x,fs=1,norm = FALSE, pl = TRUE,omit=-(1:5)){
  # Pad x with zeroes to make it's length a power of 2, i.e. length should be 2^something
  # This makes the fft faster
  oldx <- x # keep for later
  if(norm == TRUE) x <- x - mean(x)
  nfft <- (2^ceiling(log2(length(x))))
  x[(length(x)+1): nfft] <- 0
  fftx <- fft(x)
  # It's symmetrical so throw away second half. Only first 1 + nfft points are unique
  NoUnPo <- ceiling((nfft+1)/2) # Number of unique points
  fftx <- fftx[1:NoUnPo]
  # First element is DC component, last is the Nyquist component
  
  # Take magnitude of fft of x and scale the fft so that it is not a function of length
  mx <- abs(fftx) / length(x)
  # Take square of magnitude
  mx <- mx^2
  
  # As we dropped the first half of fft so multiply by 2 to keep same energy
  # DC component (first element) and Nyquist component (last element if even
  # number points (should be)) should not be multiplied by 2
  mx[2:(length(mx)-1)] <- mx[2:(length(mx)-1)] * 2
  
  f <- seq(from =0, to = NoUnPo-1) * fs/nfft # frequency axis for plot

  return(as.data.frame(list(mx = mx, f = f)))}

setMethod('spectra', signature(object='FLQuant'),function(object,...){
  adply(object,c(1,3:6),spectraFn)})
  
setMethod('spectra', signature(object='FLStock'),function(object,...){
  rbind(
  cbind(qname="stock",   spectraFLQ(stock(object))),
  cbind(qname="ssb",     spectraFLQ(ssb(  object))),
  cbind(qname="catch",   spectraFLQ(catch(object))),
  cbind(qname="rec",     spectraFLQ(rec(  object))),
  cbind(qname="f",       spectraFLQ(fbar( object))),
  cbind(qname="harvest", spectraFLQ(catch(object)/stock(object))))})
  

def blockAverage(datastream, maxBlockSize):
 
  # total number of observations in datastream
  Nobs         = len(datastream)  

# mean (expect to be "nearly" constant)
  blockMean = np.zeros(Nobs)  
# variance associated with each blockSize                   
  blockVar  = np.zeros(Nobs)        
             
  blockCtr  = 1
  
#  blockSize is number of observations/block
#  run them through all the possibilities

  for blockSize in range(1,maxBlockSize+1):

# total number of such blocks in datastream
    Nblock    = int(math.floor(Nobs/blockSize))     
    
    #print "Nblock = ", Nblock
    
    # container for parcelling block 
    obsProp   = np.zeros(Nblock)                   
 				
	# Loop to chop datastream into blocks
	# and take average
    for i in range(0,Nblock):
      ibeg = i*blockSize + 1
      iend =  ibeg + blockSize - 1
      obsProp[i-1] = np.mean(datastream[(ibeg-1):iend])
      #print "ibeg = ", ibeg              
      #print "iend = ", iend              
      #print "obsProp[i-1] = ", obsProp[i-1]              


    blockMean[blockCtr-1] = mean(obsProp)
    blockVar[blockCtr-1] = statistics.variance((obsProp)/(Nblock - 1))
    blockCtr = blockCtr + 1
    
  v=range(minBlockSize,(maxBlockSize+1))
  
  return v,blockVar,blockMean
  
  
# NOW WE CALL THE FUNCTION blockAverage()

datastream = theta[10:80]
maxBlockSize = 4
v,blockVar,blockMean = blockAverage(datastream, maxBlockSize)
print "v = ", v
print "blockVar = ", blockVar
print "blockMean = ", blockMean
print len(blockVar), len(blockMean)  , len(v)

# NOW WE PLOT OUR RESULTS --> BUT WE STILL HAVE TO CALCULATE THE FINAL ERROR



def blockAverage2(datastream, Nblock):
    
    Nobs = len(datastream) 
    maxBlockSize = int(math.floor(Nobs/Nblock)) 
    remainder = Nobs%Nblock
    
    #IF Nobs/Nblock IS NOT AN INTEGER (look for function that returns division remainder --> if renainder is not equal to zero) then add 1 to the last block 
    blockMean = np.zeros(Nobs)  
    
    for i in range(0,Nblock +1):
    
        ibeg = i*maxBlockSize
        iend = (i+1)*maxBlockSize
        blockMean = mean(datastream[ibeg:iend])
        
        
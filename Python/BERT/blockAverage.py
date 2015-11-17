def blockAverage(datastream, maxBlockSize):
 
  Nobs         = len(datastream)           # total number of observations in datastream
  minBlockSize = 1                            # min: 1 observation/block

  #print "Nobs = ", Nobs
    
  NumBlocks = maxBlockSize - minBlockSize + 1    # total number of block sizes

  blockMean = np.zeros(NumBlocks)                   # mean (expect to be "nearly" constant)
  blockVar  = np.zeros(NumBlocks)                   # variance associated with each blockSize
  blockCtr  = 1
				#
				#  blockSize is # observations/block
				#  run them through all the possibilities
				#
  #print "NumBlocks = ", NumBlocks

  for blockSize in range(minBlockSize,maxBlockSize+1):

    Nblock    = int(math.floor(Nobs/blockSize))     # total number of such blocks in datastream
    
    #print "Nblock = ", Nblock
    
    obsProp   = np.zeros(Nblock)                    # container for parcelling block 
 				
	# Loop to chop datastream into blocks
	# and take average
    for i in range(1,Nblock+1):
      ibeg = (i-1)*blockSize + 1
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
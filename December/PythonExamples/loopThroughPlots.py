fig, ax = plt.subplots(10, 10)


#Definition: 
plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, **fig_kw)
#Create a figure with a set of subplots already made.

#This utility wrapper makes it convenient to create common layouts of
#subplots, including the enclosing figure object, in a single call.

for i in range(100):
    ax = plt.subplot(10,10,i)
    ax.plot(...)
    
#####################################
#Another example:

# not tested
import math
import matplotlib.pylab as plt

Nrows = math.ceil(len(subsl) / 2.)
for i in range(len(subsl)):
    subsm = subsl[i]
    H7, subsm = sumsubdesc2(table, subsm) 
    plt.subplot(Nrows, 2, i+1)

    # do some plotting

    plt.title('Rolling 4q mean %s'%(subsm))
    
    
#####################################
# My version:
import math
import matplotlib.pylab as plt

# For 10 plots in 5 rows and 2 columns:

fig, ax = plt.subplots(5, 2)

Nrows = 5
for i in range(10):
    ax = plt.subplot(Nrows, 2, i+1)
    ax.plot(...)
    # do some plotting

    plt.title('Rolling 4q mean %s'%(subsm))
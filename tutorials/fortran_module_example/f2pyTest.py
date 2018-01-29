from pylab import *
import f2pyTest

a = 4
testmat = arange(10)

testmat_out,out2 = f2pyTest.foo(a,testmat)

print testmat_out
print out2

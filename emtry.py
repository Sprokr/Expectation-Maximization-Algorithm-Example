import scipy
import random
import Tkinter
import numpy
import math
import time
import pylab

def ellipse(ra,rb,ang,x0,y0,Nb):
	'''ra - major axis length
	rb - minor axis length
	ang - angle
	x0,y0 - position of centre of ellipse
	Nb - No. of points that make an ellipse
	'''
	xpos,ypos=x0,y0
	radm,radn=ra,rb
	an=ang
	
	co,si=scipy.cos(an),scipy.sin(an)
	the=scipy.linspace(0,2*math.pi,Nb)
	X=radm*scipy.cos(the)*co-si*radn*scipy.sin(the)+xpos
	Y=radm*scipy.cos(the)*si+co*radn*scipy.sin(the)+ypos
	return X,Y


def myfunc(mat_cordinate,mat_invariance,mat_mean):
	mat_mean=numpy.matrix.reshape(mat_mean,(2,1))
	diffmean=mat_cordinate-mat_mean
	
	det=numpy.linalg.det(mat_invariance)
	if(abs(det) < 0.0001):
		
		quit("stopped")
	
	exppower = diffmean.T * numpy.linalg.inv(mat_invariance) * diffmean	
	exppower=exppower*(-1/2)
	numerator=math.exp(exppower[0,0])
	rootdet=math.sqrt(det)
	denom=math.pi*2*rootdet
	probability=numerator/denom
	return (probability)








x1=1
x2=100
y1=1
y2=100


num_pixels=int(raw_input("Enter Number of points :- "))
num_clusters=int(raw_input("Enter Number of Clusters :- "))

points = numpy.matrix(numpy.zeros(shape=(2,num_pixels)))
mean_difference = numpy.matrix(numpy.zeros(shape=(2,num_pixels)))
mean_difference_transpose = numpy.matrix(numpy.zeros(shape=(num_pixels,2)))
covariance = numpy.matrix(numpy.zeros(shape=(2,2)))
covariance_3d = numpy.zeros(shape=(num_clusters,2,2))
cluster = numpy.matrix(numpy.zeros(shape=(1,num_pixels)))
weight = numpy.matrix(numpy.zeros(shape=(num_pixels,num_clusters)))
sum_arr=numpy.zeros(shape=(num_clusters))
avg=numpy.zeros(shape=(num_clusters))
Nk=numpy.zeros(shape=(num_clusters))
alpha_k=numpy.zeros(shape=(num_clusters))
cluster_means = numpy.matrix(numpy.zeros(shape=(2,num_clusters)))
result = numpy.matrix(numpy.zeros(shape=(num_pixels,num_clusters)))
pointsincluster = numpy.zeros(shape=(num_clusters))
weight = numpy.matrix(numpy.zeros(shape=(num_pixels,num_clusters)))
tmp = numpy.matrix(numpy.zeros(shape=(2,1)))
tmp_covar = numpy.matrix(numpy.zeros(shape=(2,2)))
prob = numpy.zeros(shape=(num_clusters))
log_likelyhood = numpy.zeros(shape=(2))
color = []
major_length = numpy.matrix(numpy.zeros(shape=(num_clusters)))
minor_length = numpy.matrix(numpy.zeros(shape=(num_clusters)))
angle = numpy.matrix(numpy.zeros(shape=(num_clusters)))
new_angle = numpy.matrix(numpy.zeros(shape=(num_clusters)))
A = numpy.matrix(numpy.zeros(shape=(2, 2)))
evecs = numpy.matrix(numpy.zeros(shape=(2, 2)))



if(num_pixels < num_clusters):
	quit("Number of pixels to be clustered can't be less then Number Of Clusters!!")
	
#Points Generator
for i in xrange(0,2):
	for j in xrange(0,num_pixels):
		points[i,j]=random.randrange(x1+1,x2)


#Assign cluster to each point randomly
i=0
for num in xrange(0,num_pixels):
	for k in xrange(0,num_clusters):
		if(i<num_pixels):
			cluster[0,i]=k+1
			i=i+1
print cluster		

for i in xrange(0,2):
	sum_arr[i]=points[i,:].sum()
	avg[i]=float(sum_arr[i])/num_pixels
	for j in xrange(0,num_pixels):
		mean_difference[i,j]=points[i,j]-avg[i]
		



mean_difference_transpose=mean_difference.T
covariance = mean_difference*mean_difference_transpose

covariance=covariance/num_pixels

for i in xrange(0,num_clusters):
	for j in xrange(0,2):
		for k in xrange(0,2):
			covariance_3d[i,j,k]=covariance[j,k]


# initialize cluster means
for i in xrange(0,2):
	for j in xrange(0,num_clusters):
		cluster_means[i,j]=random.randrange(x1+1,x2)

# create an matrix which keeps record of cluster and point in it
	for i in xrange(0,num_pixels):
		for j in xrange(1,num_clusters+1):
			if(cluster[0,i]==j):
				weight[i,j-1]=1
print weight	
for k in xrange(0,num_clusters):
	pointsincluster[k]=weight[:,k].sum()
	

count = 500
totaliteration = 0
#function Itereation
while(count > 1):
	
	# calculate alpha_k
	for i in xrange(0,num_clusters):
		alpha_k[i]=pointsincluster[i]/num_pixels
		
	
	for i in xrange (0,num_pixels):
		for j in xrange (0,num_clusters):
			result[i,j]=myfunc(points[:,i],covariance_3d[j,:,:],cluster_means[:,j] )* alpha_k[j]
		summation=result[i,:].sum()
		weight[i,:]=result[i,:]/summation
	
	for k in xrange(0, num_clusters):
		alpha_k[k] = pointsincluster[k]/num_pixels
	
	#find new cluster means	
	for k in xrange(0,num_clusters):
		tmp[0,0]=0
		tmp[1,0]=0
		for i in xrange(0,num_pixels):
			tmp[:,0]=tmp[:,0]+weight[i,k]*points[:,i]
		cluster_means[:,k]=(tmp[:,0]/pointsincluster[k])
	
	for k in xrange(0,num_clusters):
		tmp_covar[:,:]=0
		for i in xrange(0, num_pixels):
			tmp_covar[:,:] = tmp_covar[:,:] + (weight[i,k]*((points[:,i]-cluster_means[:,k])*(points[:,i]-cluster_means[:,k]).T))
		covariance_3d[k,:,:]=(tmp_covar[:,:]/pointsincluster[k])
	
	for i in xrange (0,num_pixels):
		for j in xrange (0,num_clusters):
			result[i,j]=myfunc(points[:,i],covariance_3d[j,:,:],cluster_means[:,j] )* alpha_k[j]
		summation=result[i,:].sum()
		prob[j]=numpy.log(summation)
	summation=prob.sum()
	log_likelyhood[1]=log_likelyhood[0]
	log_likelyhood[0]=summation
	diff=abs(log_likelyhood[0]-log_likelyhood[1])
	if(diff < 0.00001):
		count=1
		print "Converged"
	
	
	max_value_index = weight.argmax(axis=1)
	
	
	for i in xrange(0,num_pixels):
		#print cluster
		pos=cluster[0,i]
		cluster[0,i]=max_value_index[i,0]
		
		
	for k in xrange(0,num_clusters):
			pointsincluster[k]=weight[:,k].sum()

	count=count-1
	totaliteration=totaliteration+1
print totaliteration
	



for k in xrange(0,num_clusters):
	tmpcolor = hex(random.randint(0, 16777215))
	tmpcolor=str(tmpcolor)
	tmpcolor='#'+tmpcolor[2:8]
	color.insert(k,tmpcolor)



for k in xrange(0,num_clusters):
	A = covariance_3d[k,:,:]
	evals, evecs = numpy.linalg.eig(A)
	evecs_det = numpy.linalg.det(evecs)
	if(evals[0] >= evals[1]):
		major_length[0,k]=2*math.sqrt(evals[0])
		minor_length[0,k]=2*math.sqrt(evals[1])
		angle[0,k] = math.degrees(math.atan(float(evecs[1,0]/evecs[0,0])))
	else:
		major_length[0,k]=2*math.sqrt(evals[1])
		minor_length[0,k]=2*math.sqrt(evals[0])
		angle[0,k] = math.degrees(math.atan(float(evecs[1,1]/evecs[0,1])))
	if(float(evecs[1,1]/evecs[0,0]) >= 0):		
		new_angle[0,k] = angle[0,k]
	else:
		new_angle[0,k] = 90-angle[0,k]
	X,Y=ellipse(major_length[0,k],minor_length[0,k],new_angle[0,k]*(math.pi/180),cluster_means[0,k],cluster_means[1,k],100)
	pylab.plot(X,Y,color[k-1],ms=1) # blue ellipse

		

for k in xrange(0,num_clusters):
	cir = pylab.Rectangle((cluster_means[0,k],cluster_means[1,k]), 3,3,  fc=color[k-1])
	pylab.gca().add_patch(cir)
	
for i in xrange(0,num_pixels):
	x=points[0,i]
	y=points[1,i]
	k = int(cluster[0,i])-1
	color[k]=color[k]
	cir = pylab.Circle((x,y), radius=1,  fc=color[k])
	pylab.gca().add_patch(cir)
	

pylab.show()
pylab.show()

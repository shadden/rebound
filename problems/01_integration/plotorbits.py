data = loadtxt("orbits.txt")
p1dat = data[::2]
p2dat = data[1::2]
plot( p1dat[:,0],p1dat[:,1]*(1+p1dat[:,2]),'k-' )
plot( p1dat[:,0],p1dat[:,1]*(1-p1dat[:,2]) ,'k-')
plot( p1dat[:,0],p1dat[:,1] ,'k-')
plot( p2dat[:,0],p2dat[:,1]*(1+p2dat[:,2]) ,'k-')
plot( p2dat[:,0],p2dat[:,1]*(1-p2dat[:,2]) ,'k-')
plot( p2dat[:,0],p2dat[:,1] ,'k-')
show()

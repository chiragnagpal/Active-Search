
# coding: utf-8

# In[3]:

class AS:
        
    
    def __init__(self, A, N,eta, w, pi, lab_pos, lab_neg,lab_warm,amg):
        
         
        import numpy as np
        from scipy.sparse import coo_matrix       

        
        self.lab_pos = lab_pos
        
        self.lab_neg = lab_neg
        
        self.lab_warm = lab_warm
        
        self.amg = amg
 
        self.N = N
        
        self.A = A

#       unlab = list( ( (set(range(N)) - set(lab_pos)) -set(lab_neg)) - set(lab_ign)  )

        lab = set(lab_pos) | set(lab_neg)| set(lab_warm)

        unlab = set(range(N)) - lab

        #convert A to matrix
        
        self.itr = -1
          
       # setting D

        self.D = np.asmatrix(np.diag(np.squeeze(np.asarray(self.A.sum(axis =1)))))


        #print "D is set"
        #D is set

        #setting B
        self.B = np.zeros(N)
        for i in range(N):

            if i in lab:

                self.B[i] = eta/(1-eta)
            if i in unlab:
                self.B[i] = 1./(1+w)

                
        #self.B = np.asmatrix(self.B)
        #print "B is set"        
        #B is set
        #setting I
        self.I = np.identity(N)
        #I is set
        #setting y
        
        self.y = np.zeros(N)
        
        for i in range(N):
            
            if i in lab_warm:

                self.y[i] = 0.5
            
            if i in lab_pos:

                self.y[i] = 1

            if i in lab_neg:

                self.y[i] = 0

            if i in unlab:
                self.y[i] = pi

	#self.y = np.asmatrix(self.y)
	
        #y is set

        #print "y is set"
    
	#convert to CSR
        print "converting to CSR"
        
        
        self.eta = eta
        self.w = w
        self.pi = pi
        #self.y = np.asmatrix(self.y)
        if True:
        
        #if  self.amg:

            self.A = coo_matrix(self.A)
            self.D = coo_matrix(self.D)
            self.I = coo_matrix(np.identity(N))
#           self.B = coo_matrix(self.B)
#           self.y = coo_matrix(self.y)

            self.A = self.A.tocsr()
#           self.B = self.B.tocsr()
            self.D = self.D.tocsr()
#           self.y = self.y.tocsr()
            self.I = self.I.tocsr()

        print 'A:', self.A.shape
        print 'B:', self.B.shape
        print 'D:', self.D.shape
        print 'y:', self.y.shape
        print 'I:', self.I.shape

        print "converted"
    
    
    def update(self, index, label):
        
        #self.y = self.y.todense()
        
        #self.B = self.B.todense()
        
        
        if label == 'pos':
            
            self.lab_pos = self.lab_pos+[index]
            
            self.y[index] = 1

        
        else:
            
            self.lab_neg = self.lab_neg+[index,]
            
            self.y[index] = 0

        
        from scipy.sparse import csr_matrix    
        
        
        
        self.B[index] = self.eta/(1.-self.eta)
            
        
        self.setCH()
         
    def search(self, itr=False):
        import numpy as np
    

        
        if self.amg:
            
            if self.itr<0:
                
                self.itr+=1
            
                from pyamg import ruge_stuben_solver
                from scipy.sparse import diags

                tempB = diags([self.B.tolist()] , [0])
                tempy = self.y

                #self.C =  (self.D *  (self.I + self.B) ) -self.A

                self.C =  (self.D *  (self.I + tempB) ) -self.A


                #self.H = self.D *  (tempB * self.y.T)

                self.H = self.D * (tempB * tempy)

                #self.H = self.H.todense()

                self.H = np.asarray(self.H)


                ml = ruge_stuben_solver(self.C)

                self.ml = ml

                if  itr:
                    solution = ml.solve(self.H, tol=0.5e-4, return_residuals=True)
                    f = np.squeeze(solution[0])
                    x = solution[1]

                    return (f, x)

                else:
                    f = np.squeeze(ml.solve(self.H, tol=0.5e-4))
                    self.out = f
                    return f


            
            else:
                
                self.ml = self.create_solver()
                
                self.out = self.ml.solve(self.H, tol=0.5e-4,x0=self.out,return_residuals=False)
                
                return self.out

        else:

            from scipy.sparse.linalg import spsolve

            self.C =  (self.D *  (self.I + self.B) ) -self.A

            self.H = self.D *  np.dot(np.diag(self.B) , self.y.T)


            f = np.squeeze(spsolve(self.C, self.H))

            return f             

    
    def setCH(self):
        
        import numpy as np
        from pyamg import ruge_stuben_solver
        from scipy.sparse import diags
        
        tempB = diags([self.B.tolist()] , [0])
        tempy = self.y
            
        self.C =  (self.D *  (self.I + tempB) ) -self.A

            
            
        self.H = self.D * (tempB * tempy)

            
        #self.H = self.H.todense()
            
        self.H = np.asarray(self.H)
        
    def create_solver(self):
        
        A = self.C
        
        ml = self.ml

        ml.levels[0].A = A
    
        for i in range(len(ml.levels)-1):
        
            ml.levels[i+1].A = ml.levels[i].R * ml.levels[i].A * ml.levels[i].P 
        
        return ml


# In[4]:

def largest_index(f, l):

    l = set(l)
    
    
    temp = -1.1
    idx = -1

    for i in range(len(f)):

        if i not in l:

            if f[i]>temp:
                
                temp = f[i]
                idx = i
            
    return idx        


# In[ ]:




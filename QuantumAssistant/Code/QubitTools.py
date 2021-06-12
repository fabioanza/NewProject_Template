import numpy as np

class Qubit:
    def __init__(self,basis):
        #Initialized with a pair of orthogonal kets.
        v0,v1 = np.array(basis[0],dtype=complex),np.array(basis[1],dtype=complex)
        if np.isclose(np.vdot(v0,v1),0,atol=10**(-8))==True:
            self.basis=basis
        else:
            print('The vectors provided do not form a basis. Using computational basis')
            self.basis=np.eye(2)


###############################################################################
###################### GENERIC FUNCTIONS ALWAYS USEFUL ########################
###############################################################################
    def pauli_matrices(self,which):
        if which.lower()=='x':
            return np.array([[0,1],[1,0]],dtype=complex)
        elif which.lower()=='y':
            return np.array([[0,-1j],[1j,0]],dtype=complex)
        elif which.lower()=='z':
            return np.array([[1,0],[0,-1]],dtype=complex)
        else:
            print('Quitting. Argument {val} can only be `x`, `y` or `z`'.format(val=which))

    def scalar_product(self,psi1,psi2):
        return np.vdot(psi1,psi2)

    def regularize_ket(self,psi):
        #normalization and first phase equal to zero.
        return psi*np.exp(-1j*(np.angle(psi[0])))/np.sqrt(np.real(self.scalar_product(psi,psi)))

    def pauli_kets(self,which,excitation):
        AUX = self.pauli_matrices(which)
        e_val, e_vec = np.linalg.eigh(AUX)
        if excitation==0:
            return self.regularize_ket(e_vec[0])
        elif excitation==1:
            return self.regularize_ket(e_vec[1])
        else:
            print('Quitting. Excitation {val} can only be 0 or 1'.format(val=excitation))

    def pauli_vector(self,v0,vx,vy,vz):
        return v0*np.eye(2)+vx*self.pauli_matrices('x')+vy*self.pauli_matrices('y')+vz*self.pauli_matrices('z')

    def build_hamiltonian(self,h0,hx,hy,hz):
        self.hamiltonian = self.pauli_vector(h0,hx,hy,hz)

    def isbasis(self,vec0,vec1):
        #Boolean output to check if two vectors are a basis.
        v0,v1 = self.regularize_ket(vec0),self.regularize_ket(vec1)
        assert np.isclose(self.scalar_product(v0,v0),1,atol=10**(-8))==True, 'First argument is not a normalizable state'
        assert np.isclose(self.scalar_product(v1,v1),1,atol=10**(-8))==True, 'Second argument is not a normalizable state'
        if np.isclose(self.scalar_product(v0,v1),0,atol=10**(-8))==True:
            return True
        else:
            return False

    class State:
        def __init__(self,which,*args):
            if which=='ket':
                self.ket = np.array(args)
            if which=='p_phi':
                self.ket =

    def psi_to_pphi(self,psi,vec0,vec1):
        #Extract the (p,phi) representation of a generic ket |psi> in the basis vec0,vec1
        v0,v1 = self.regularize_ket(vec0), self.regularize_ket(vec1)
        assert self.isbasis(vec0,vec1)==True, 'Second and Third arguments do not form a basis. Quitting.'
        assert len(psi)==2, 'Length {val} should be 2.'.format(val=len(psi))
        return (np.abs(self.scalar_product(psi,vec1))**2,np.angle(self.scalar_product(psi,vec1))-np.angle(self.scalar_product(psi,vec0)))


    def pphi_to_psi(pphi,vec0,vec1):
        #Extract the ket representation |psi(p,phi)> of a CP1 point (p,phi)
        assert self.isbasis(vec0,vec1)==True, 'Second and Third arguments do not form a basis. Quitting.'
        if np.isclose(pphi[0]*(1-pphi[0]),0,atol=10**(-8))==True:
            return np.sqrt(1-pphi[0])*vec0+np.sqrt(pphi[0])*vec1
        else:
            return np.sqrt(1-pphi[0])*vec0+np.sqrt(pphi[0])*np.exp(1j*pphi[1])*vec1

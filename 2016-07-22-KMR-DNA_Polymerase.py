#get_ipython().magic(u'matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint
from sys import exit

# In[2]:

def plot_template(template,primer=0,yt=1):
    """
    plot template using a number line
    """
    plt_min = template[0] - 2
    plt_max = template[-1] + 2
    plt.xlim((plt_min,plt_max))
    plt.ylim((yt-1,2))
    if primer==1:
        plt_min = plt_min + 2 -0.2
    x = np.arange(plt_min,plt_max,1)
    y = yt*np.ones(len(x))
    plt.plot(x,y,'-k',lw=2)
    y2 = yt*np.ones(len(template))
    plt.plot(template,y2,'.k',ms=30)
    plt.axis('off')
    


# In[3]:

def reaction1(template, primer):
    """
    Megaprimer synthesis
    """
    temp = template[template > primer[0]]
    return np.hstack((primer,temp))


def reaction2(template, primer):
    """
    Mutagenesis PCR
    """
    temp = template[template<primer[0]]
    return np.hstack((temp,primer))


def set_primers(Nmut=4):
    """
    Set up initial mixture elements
    """
    primers = []
    for i in range(0,Nmut):
        ele = [i]
        primers.append(ele) 
    primers = np.array(primers)
    return primers



# In[7]:


def do_round(mixture, primers, Nmut):
    """
    Add initial_mixture elements one-by-one
    to the current mixture and get new products
    """

    total = np.zeros(Nmut)
    out_mixture = []
    for template in mixture:
        for primer in primers:
            product1 = reaction1(template,primer)
            for template2 in mixture:
                product2 = reaction2(template2,product1)
                out_mixture.append(product2)
                total[len(product2)-1] += 1
    
    Ntotal = np.sum(total)
    total = total/Ntotal
    print Ntotal, total
    return total, Ntotal, out_mixture


def test_reaction1():
    """
    function to test reaction1
    with random primers and templates
    """
    Nmut = randint(4,10,1)
    a = np.arange(0,Nmut)
    
    t = randint(0,2,Nmut)
    template = a[t==1]

    primer = randint(0,Nmut,1)
    if len(primer) == 0 or len(template)==0:
        return 0
    
    print template, "  ", primer, "  ", reaction1(template,primer)
    

def test_reaction2():
    """
    function to test reaction2
    with random primers and templates
    """
    Nmut = randint(4,10,1)

    a = np.arange(0,Nmut)
    t = randint(0,2,Nmut)
    template = a[t==1]

    a = np.arange(0,Nmut)
    p = randint(0,2,Nmut)
    primer = a[p==1]
    if len(primer) == 0 or len(template)==0:
        return 0

    print template, "  ", primer, "  ", reaction2(template,primer)
    

def do_tests():
    """
    Do reaction1, reaction2 tests
    """
    print "Reaction1"
    print "Template  ", "Primer  ", "Product"    
    for i in range(0,30):
        test_reaction1()

    print "Reaction2"
    print "Template  ", "Primer  ", "Product"    
    for i in range(0,30):
        test_reaction2()



Nmut = 8
nrounds = 3

primers = set_primers(Nmut)
initial_mixture = np.copy(primers)

for i in np.arange(0,nrounds-1):
    total, Ntotal, next_mixture = do_round(initial_mixture, primers, Nmut)
    initial_mixture = list(next_mixture)







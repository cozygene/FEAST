import numpy as np
from scipy.optimize import nnls
#from scipy.optimize import minimize
from random import sample
import matplotlib.pyplot as plt
from collections import Counter
from operator import itemgetter
import timeit
PATHTODATA="~/FEAST_Data/"
DATAFILE="M3"
NSOURCES=20 ###CHANGE ACCOMPANYING n IN unknown_initialize function usually 60% of NSOURCES
MIXINGITERATIONS=50
SAMPLEITERATIONS=10
EMITERATIONS=100
INCLUDEEPSILON=True
CLSBOOTSTRAP=True
EPSILON=.01
COVERAGE=10000
#######################



'''
    We would like to acknowledge Knight et al. for the jsd and rarefy functions.
    This code is from Knight et al., it is no way original, we just needed
    an implementation in python.
'''
def jsdmatrix(X):
    d = np.zeros(X.shape)
    for i in np.arange(d.shape[0]-1):
        for j in np.arange(i+1, d.shape[0]):
            d[i, j] = jsd(X[i, :], X[j, :])
            d[j, i] = d[i, j]
    return d

def jsd(p, q):
    m = (p + q) / 2.
    return (kld(p, m) + kld(q, m))/2.

def kld(p, q):
    nonzero = np.logical_and(p > 0, q > 0)
    return np.sum(p[nonzero] * np.log2(p[nonzero]/q[nonzero]))


def rarefy(X, newdepth):
    nrow = (X.shape)[0]
    ncol = (X.shape)[1]
    for i in range(nrow):
        if np.sum(X[i,:]) > newdepth:
            s = np.random.choice(np.arange(1, 1+ncol), size=newdepth, replace=True, p=X[i, :]/np.sum(X[i, :]))
            X[i, :] = np.histogram(s, bins=np.arange(0.5, ncol+1+.5, 1))[0]
    return X

def iscor(vec, n=3):
    if np.sum(vec > n):
        return 1
    return 0

def unknown_initialize(sources, sink, n_sources):
    ###if 60% of the sources have it, call it cor
    n = np.round(.6*NSOURCES)
    unknown_source = np.zeros(sink.shape[0])
    sum_sources = np.sum(sources, axis=0)
    sinksumdiff = np.subtract(sink,sum_sources)
    negmask = sinksumdiff < 0.0
    sinksumdiff[negmask] = 0.0
    cor_indmask = sources > 0.0
    cor_indmasksum = np.sum(cor_indmask, axis=0)
    cor_ind = cor_indmasksum > n
    otu_idx = np.asarray(range(cor_ind.shape[0]))
    cor_index = otu_idx[cor_ind]
    cor_abundance = np.round(np.apply_along_axis(np.min, 0, sources[:, cor_index]))
    unknown_source[cor_index] = cor_abundance
    return unknown_source

def minfn(x, A, b):
    return np.linalg.norm(A.dot(x) - b)

def getR2(y, fx):
    ybar = np.sum(y) / y.shape[0]
    #ssreg = np.sum((fx - ybar) ** 2)
    ssres = float(np.sum(y - fx))**2
    sstot = float(np.sum(y - ybar))**2
    return 1.0 - ssres / sstot

def create_alphas(num_sources, n):
    m= []
    if n == 1:
        index = sample(range(num_sources), 1)
        m_1 = np.random.uniform(.6, .9)
        resid = 1 - m_1
        other_ms = resid/(num_sources-1)
        m = np.nan + np.zeros(num_sources)
        m[index] = m_1
        m[np.isnan(m)] = other_ms
        #return(m)
    elif n == 2:
        index = sample(range(num_sources), 2)
        m_1 = np.random.uniform(.1, .2)
        m_2 = np.random.uniform(.4, .5)
        resid = 1 - m_1 - m_2
        other_ms = resid/(num_sources-2)
        m = np.nan + np.zeros(num_sources)
        m[index] = [m_1, m_2]
        m[np.isnan(m)] = other_ms
        #return(m)
    elif n == 3:
        index = sample(range(num_sources), 3)
        m_1 = np.random.uniform(.1, .5)
        m_2 = np.random.uniform(.2, .25)
        m_3 = np.random.uniform(.1, .15)
        resid = 1 - m_1 - m_2 - m_3
        other_ms = resid/(num_sources-3)
        m = np.nan + np.zeros(num_sources)
        m[index] = [m_1, m_2, m_3]
        m[np.isnan(m)] = other_ms
        #return(m)
    #epprop = EPSILON#abs(np.random.uniform(EPSILON-.15, EPSILON+.15))
    #epval = epprop*np.sum(m)/(1-epprop)
    m = np.asarray(m)
    epval = EPSILON*np.sum(m)/(1-EPSILON)
    m = np.append(m, epval)
    '''
    subsum = 0
    while subsum < abs(EPSILON - 1e-03):
        tosub = EPSILON - subsum
        tosub = float(tosub)/float(num_sources+1)
        mask = m > tosub
        m[mask] -= tosub
        subsum += m[mask].shape[0]*tosub
    #m /= np.sum(m)
    m = np.append(m, EPSILON)
    '''
    return m/np.sum(m)
    #return(np.zeros(num_sources) + 1.0/num_sources)

def h(x):
    y=x[x>0]
    return(-sum(y*np.log(y)))

def JSD(p,q):
    return(h(np.dot(q,p)) - np.dot(q, np.apply_along_axis(h, 1, p)))

def M_noise(sources, alphas, sink, observed):
    rel_sink = sink / float(np.sum(sink))
    if np.sum(sources[0]) > 0:
        rowsums = np.sum(sources, axis=1)
        for i in range(sources.shape[0]):
            sources[i] /= float(rowsums[i])
    newgammas = np.multiply(alphas.reshape(alphas.shape[0], 1), sink) + observed
    newgammas /= np.sum(newgammas, axis=0)
    mask = np.isnan(newgammas)
    newgammas[mask] = 0.0
    alpha_nums = np.multiply(sink, newgammas)
    alpha_nums = np.multiply(alphas.reshape(alphas.shape[0], 1), alpha_nums)
    alpha_nums = np.sum(alpha_nums, axis=1)
    alpha_denom = np.multiply(alphas.reshape(alphas.shape[0], 1), newgammas)
    alpha_denom = np.sum(alpha_denom)
    newalphas = np.divide(alpha_nums, alpha_denom)
    newalphas = np.divide(newalphas, np.sum(newalphas))
    return [newgammas, newalphas]

def E(sources, alphas):
    nums = np.multiply(alphas.reshape(alphas.shape[0], 1), sources)
    nums = np.sum(nums, axis=1)
    denom = np.sum(nums)
    return nums / float(denom)

def do_EM_noise(sources, alphas, sink, em_itr, observed):
    curalphas = alphas; newalphas = alphas;
    m_guesses = [alphas[0]]
    for i in range(em_itr):
        curalphas = E(sources, newalphas)
        tmp = M_noise(sources, curalphas, sink, observed)
        newalphas = tmp[1]
        sources = tmp[0]
        m_guesses.append(newalphas[0])
        if abs(m_guesses[len(m_guesses)-1]-m_guesses[len(m_guesses)-2]<=10**-6):
            break
    return [sources, newalphas]

def M_basic(sources, alphas, sink):
    XOs = np.multiply(sources, sink)
    AOs = np.multiply(alphas.reshape(alphas.shape[0], 1), sources)
    newAs = np.empty(alphas.shape[0])
    for i in range(alphas.shape[0]):
        newA = np.multiply(alphas.reshape(alphas.shape[0], 1)[i],
                           np.divide(XOs[i], np.sum(AOs, axis=0)))
        mask = np.isnan(newA)
        newA[mask] = 0.0
        newA = np.sum(newA, axis=0)
        newAs[i] = newA
    return(np.divide(newAs, sum(newAs)))

def do_EM_basic(sources, alphas, sink, em_itr):
    curalphas = alphas; newalphas = alphas;
    m_guesses = [alphas[0]]
    for itr in range(em_itr):
        curalphas = E(sources, newalphas)
        newalphas = M_basic(sources, curalphas, sink)
        m_guesses.append(newalphas[0])
        if abs(m_guesses[len(m_guesses)-1] - m_guesses[len(m_guesses) - 2] <= 10 ** -6):
            break
    return newalphas


def run_sim(dat, itr, eps, bs, K):
    loadstr = PATHTODATA+dat+"_stripped.txt"
    source = np.loadtxt(loadstr, dtype="float")
    cls_R2_values = []; cls_runtimes = [];
    em_basic_R2_values = []; em_basic_runtimes = [];
    em_noise_R2_values = []; em_noise_runtimes = [];
    JS_values = []
    '''
    epsilon=0
    alpha_unobs=0.0
    if eps:
        K += 1
        alpha_unobs = np.random.uniform(.045, .055)
    '''
    for it in range(itr):
        #print "iteration: " + str(it)
        gammamask = np.random.randint(0, source.shape[0], K)
        gamma = source[gammamask]
        gamma = rarefy(gamma, COVERAGE)
        unksource = []
        if INCLUDEEPSILON:
            unksource = source[np.random.randint(0, source.shape[0], 1)]
        unksource = rarefy(unksource, COVERAGE)
        C = np.sum(gamma, axis=1)[1]
        #print C
        OTU_count = gamma.shape[1]
        y = np.empty((K, OTU_count))
        gamma /= float(C)
        for k in range(K):
            y[k] = (np.random.multinomial(C, gamma[k], size=1))


        alphas = np.empty((MIXINGITERATIONS, K+1))
        x = np.empty((MIXINGITERATIONS, OTU_count))
        beta = np.empty((MIXINGITERATIONS, OTU_count))
        for m in range(MIXINGITERATIONS):
            alphas[m] = create_alphas(K, n=3)
            beta[m] = np.dot(alphas[[m]], np.vstack([gamma, unksource/float(C)]))
            x[m] = np.random.multinomial(C, beta[m], size=1)  ##this is the observed sink
            #x[m] = beta[m]*float(C)
            #x[m] = np.random.poisson(lam=x[m])

        #weights = np.zeros(K) + 1.0/ float(K)
        #JS_values.append(JSD(gamma, weights))
        JSDMatrix = jsdmatrix(gamma[np.arange(NSOURCES), :])
        #JSDMatrix /= COVERAGE
        JS = np.mean(JSDMatrix[np.nonzero(JSDMatrix)])
        #print JS
        JS_values.append(JS)
        bs_y = None
        if bs:
            bs_y = np.empty(gamma.shape)
            for k in range(K):
                bs_y[k] = np.sum(np.random.multinomial(C, gamma[k], size=100),axis=0)
                bs_y[k] /= float((C*100))
        ###CLS
        #cons = {'type': 'eq', 'fun': lambda x: np.sum(x)-1}
        #bounds = ((0.0, 1.0),)*K
        cls_R2 = []; em_basic_R2 = []; em_noise_R2 = [];
        cls_rts = []; em_basic_rts = []; em_noise_rts = [];
        for m in range(MIXINGITERATIONS):
            unknowny = unknown_initialize(y, x[m], NSOURCES)
            cls_start_time = timeit.default_timer()
            A = np.vstack((y, unknowny)).T  #y.T  #np.vstack((gamma.T, np.ones(K)))
            b = x[m]   #np.append(beta[0], 1)
            A /= float(C)
            b /= float(C)
            cls_alphas = nnls(A, b)[0]
            cls_rts.append(timeit.default_timer() - cls_start_time)
            if np.any(np.greater(cls_alphas, 1.0)):
                cls_alphas = cls_alphas / float(np.sum(cls_alphas))
            cls_R2.append((np.corrcoef(alphas[m], cls_alphas)**2)[0,1])
            em_basic_start_time = timeit.default_timer()
            if INCLUDEEPSILON:
                em_basic_alphas = do_EM_basic(np.vstack((y, unknowny)), cls_alphas, x[[m]], EMITERATIONS)
            else:
                em_basic_alphas = do_EM_basic(y, cls_alphas, x[[m]], EMITERATIONS)
            em_basic_rts.append(timeit.default_timer() - em_basic_start_time)
            em_noise_start_time = timeit.default_timer()
            if INCLUDEEPSILON:
                em_noise_alphas = do_EM_noise(np.vstack((y, unknowny)), cls_alphas, x[[m]], EMITERATIONS, np.vstack((y, unknowny)))[1]
            else:
                em_noise_alphas = do_EM_noise(y, cls_alphas, x[[m]], EMITERATIONS, y)[1]
            em_noise_rts.append(timeit.default_timer() - em_noise_start_time)
            em_basic_R2.append((np.corrcoef(alphas[m], em_basic_alphas)**2)[0,1])
            em_noise_R2.append((np.corrcoef(alphas[m], em_noise_alphas)**2)[0,1])
        cls_R2 = np.asarray(cls_R2); em_basic_R2 = np.asarray(em_basic_R2); em_noise_R2 = np.asarray(em_noise_R2);
        cls_R2_values.append(np.mean(cls_R2[~np.isnan(cls_R2)]))
        cls_runtimes.append(np.mean(cls_rts))
        em_basic_R2_values.append(np.mean(em_basic_R2[~np.isnan(em_basic_R2)]))
        em_basic_runtimes.append(np.mean(em_basic_rts))
        em_noise_R2_values.append(np.mean(em_noise_R2[~np.isnan(em_noise_R2)]))
        em_noise_runtimes.append(np.mean(em_noise_rts))
    fig, ax = plt.subplots()
    c = Counter(JS_values).items()
    c.sort(key=itemgetter(0))
    ordered = zip(JS_values, cls_R2_values, em_basic_R2_values, em_noise_R2_values)
    ordered = sorted(ordered)
    JS_values = [x[0] for x in ordered]
    cls_R2_values = [x[1] for x in ordered]
    em_basic_R2_values = [x[2] for x in ordered]
    em_noise_R2_values = [x[3] for x in ordered]
    ax.plot(JS_values, cls_R2_values, label='CLS R2')
    ax.plot(JS_values, em_basic_R2_values, label='EM basic R2')
    ax.plot(JS_values, em_noise_R2_values, label='EM w/ noise R2')
    legend = ax.legend(loc='best', shadow=True)
    print [str(x) for x in cls_R2_values]
    print [str(x) for x in em_basic_R2_values]
    print [str(x) for x in em_noise_R2_values]
    if np.array_equal(cls_R2_values, em_basic_R2_values):
        print "CLS and EM basic R2 equal"
    if np.array_equal(cls_R2_values, em_noise_R2_values):
        print "CLS and EM noise R2 equal"
    if np.array_equal(cls_R2_values, em_basic_R2_values):
        print "noise and basic R2 equal"
    #plt.ylim([0.00, 1.00])
    plt.grid()
    plt.xlabel("JS value")
    plt.ylabel("R^2 value")
    #fig1, axn = plt.subplots()
    #axn.xaxis.set_visible(False)
    #axn.yaxis.set_visible(False)
    cls_runtimes = np.mean(cls_runtimes)
    em_noise_runtimes = np.mean(em_noise_runtimes)
    em_basic_runtimes = np.mean(em_basic_runtimes)
    #collabel = ("CLS", "EM basic", "EM noise")
    #plt.table(cellText=np.asarray([cls_runtimes, em_basic_runtimes, em_noise_runtimes]),
    #          colLabels=collabel, cellLoc='center', rowLoc='center', loc='bottom') #bbox=[0.25, -5, .9, .3])
    #plt.subplots_adjust(left=.2, bottom=.1)
    #axn.table(cellText=np.asarray([cls_runtimes, em_basic_runtimes, em_noise_runtimes]), colLabels=collabel,loc='top')
    #print [str(x) for x in cls_runtimes]
    #print [str(x) for x in em_basic_runtimes]
    #print [str(x) for x in em_noise_runtimes]
    print cls_runtimes
    print em_basic_runtimes
    print em_noise_runtimes
    #table = ax.table(cellText=np.asarray([cls_runtimes, em_basic_runtimes, em_noise_runtimes]), colLabels=collabel, loc='center')
    epstring = ""
    if INCLUDEEPSILON:
        epstring = "_epsilon" + str(EPSILON).split('.')[1] + "_"
    svstring = PATHTODATA + dat + "_"+str(NSOURCES) + "sources"+ "_mixing"+str(MIXINGITERATIONS) +"_sampling" + str(SAMPLEITERATIONS) +"_clsinit_"+ epstring
    plt.savefig(svstring)
    with open(svstring + "_runtimes.txt", 'w') as wfile:
        wfile.write('CLS\t')
        wfile.write("EMbasic\t")
        wfile.write("EMnoise\t")
        wfile.write('\n')
        wfile.write('{}\t'.format(cls_runtimes))
        wfile.write('{}\t'.format(em_basic_runtimes))
        wfile.write('{}\t'.format(em_noise_runtimes))
    #plt.show()
    #savestr




run_sim(DATAFILE, SAMPLEITERATIONS, INCLUDEEPSILON, CLSBOOTSTRAP, NSOURCES)
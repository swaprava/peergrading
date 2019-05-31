#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 12:17:40 2019

This code finds the counterexample that reliability monotone cannot be guaranteed for k neq 1/2. 

Part of the computer aided proof of the uniqueness result.

@author: swaprava
"""

import numpy as np
import scipy as sp
from scipy.integrate import quad

def integrand(x, k, sigma):

    temp_1 = np.sum( (m_k - np.mean(m_k))**2.0 ) / float(len(m_k))
    K = np.dot(np.sqrt(tau_l_hat),n_j - np.mean(n_k, axis=1)) - (sigma * (m_j - np.mean(m_k))) / (temp_1 * sigma**2)**k

    toReturn = x * (1/np.sqrt(2*np.pi)) * np.exp(-(x - K)**2)
    return toReturn


def func(k, sigma, htype='square'):
    
#    This will be needed for the last part of the payment part
#    Here we compute only the W_j^* part
    
    Z_minus_i = np.sqrt(gamma) * (mu-y_j) + np.dot(np.sqrt(tau_l_hat),difference_term_l)
    X_minus_i = np.sqrt(gamma) + np.sum(np.sqrt(tau_l_hat))
    temp_1 = np.sum( (m_k - np.mean(m_k))**2.0 ) / float(len(m_k))

    if htype == 'square':
        t_i_j = -( ( Z_minus_i * (temp_1 * sigma**2)**k + sigma * (m_j - np.mean(m_k)) ) / \
        ( X_minus_i * (temp_1 * sigma**2)**k + 1) )**2

    tau_l_hat_wo_one = tau_l_hat[1:]
    difference_term_l_wo_one = difference_term_l[1:]
    Z_minus_i_wo_one = np.sqrt(gamma) * (mu-y_j) + np.dot(np.sqrt(tau_l_hat_wo_one),difference_term_l_wo_one)
    X_minus_i_wo_one = np.sqrt(gamma) + np.sum(np.sqrt(tau_l_hat_wo_one))

    if htype == 'square':
        t_i_j_minus_t_1_j = -( ( Z_minus_i_wo_one * (temp_1 * sigma**2)**k + sigma * (m_j - np.mean(m_k)) ) / ( X_minus_i_wo_one * (temp_1 * sigma**2)**k + 1) )**2

    tau_l_hat_wo_last = tau_l_hat[:-1]
    difference_term_l_wo_last = difference_term_l[:-1]
    Z_minus_i_wo_last = np.sqrt(gamma) * (mu-y_j) + np.dot(np.sqrt(tau_l_hat_wo_last),difference_term_l_wo_last)
    X_minus_i_wo_last = np.sqrt(gamma) + np.sum(np.sqrt(tau_l_hat_wo_last))

    if htype == 'square':
        t_i_j_minus_t_3_j = -( ( Z_minus_i_wo_last * (temp_1 * sigma**2)**k + sigma * (m_j - np.mean(m_k)) ) / ( X_minus_i_wo_last * (temp_1 * sigma**2)**k + 1) )**2

    term_2 = w * (t_i_j_minus_t_1_j + t_i_j_minus_t_3_j)
    term_3 = (1 - w * 2) * t_i_j

    integralValue, error = quad(integrand, 0, np.inf, args=(k,sigma))
    multiplier = 1.0 / (X_minus_i + 1.0 / ((temp_1 * sigma**2)**k) )

    term_1 = multiplier * integralValue

    # print 'term_1 =', term_1, '\nterm_2 =', term_2, '\nterm_3 =', term_3
    # print type(term_1), type(term_2), type(term_3)
    
    return -w * term_1 + term_2 + term_3
    
    
if __name__ == '__main__':
    
    maxNumOfRuns = 20000
#    Generate all the variables here
    lenOfProbe = 3
    lenOfNonProbe = 3
    sigma_1, sigma_2 = 1, 2
    # kVector = [0.15, 0.25, 0.35, 0.45, 0.5, 0.55, 0.65, 0.75, 0.85, 1.0, 1.25, 1.5, 1.75, 2.0]
    kVector = np.linspace(0.47, 0.55, 17)
    maxPositiveChange = np.zeros(len(kVector))
    parametersForCounterexample = {}

    # k = 0.75
    w = 1/2.0 # denom has to be >= 2.0
    
    '''
    for run in range(maxNumOfRuns):
        m_k = np.random.normal(0,1,size=lenOfProbe) 
        m_j = np.random.randn()
        tau_l_hat = np.random.uniform(0,1,size=lenOfNonProbe-1) # for other graders
        # var_l = np.random.uniform(0,1,size=lenOfProbe)
        n_k = []
        for i in range(lenOfNonProbe-1):
            temp = []
            for j in range(lenOfProbe):
                temp.append(np.random.normal(0, np.random.rand()))
            n_k.append(temp)
            
        n_k = np.array(n_k) # rows = l indices, columns = probes for each of them
        n_j = []
        for i in range(lenOfNonProbe-1):
            n_j.append(np.random.normal(0, np.random.rand()))
        n_j = np.array(n_j)

        gamma = np.random.rand()
        mu = np.random.randn()
        difference_term_l = sp.random.normal(0,1,size=lenOfNonProbe-1)
        y_j = np.random.randn()

        stringOfParameters = 'm_k = ' +repr(m_k)+ '\nm_j = ' +repr(m_j)+ \
            '\ntau_l_hat = ' +repr(tau_l_hat)+ '\nn_k = ' +repr(n_k)+ \
            '\nn_j = ' +repr(n_j)+ '\ngamma = ' +repr(gamma)+ '\nmu = ' +repr(mu)+ \
            '\ndifference_term_l = ' +repr(difference_term_l)+ '\ny_j = ' +repr(y_j)

        for kIndex in xrange(len(kVector)):

            k = kVector[kIndex]
            change = func(k,sigma_2, 'square') - func(k,sigma_1, 'square')

            if change > 0:
                print 'counterexample found for k = ', k, 'change =', change
                if change > maxPositiveChange[kIndex]:
                    maxPositiveChange[kIndex] = change
                    parametersForCounterexample[k] = stringOfParameters


    f = open('uniqueness_results_zoomed_near_half.txt', 'w')

    print 'kVector =', kVector
    print 'maxPositiveChange =', maxPositiveChange

    f.write('kVector = ' +repr(kVector)+ '\n')
    f.write('maxPositiveChange = ' +repr(maxPositiveChange)+ '\n')

    print 'parametersForCounterexample =\n'
    f.write('\n==========================\nSummary of Counterexamples\n==========================\n')
    for kIndex in xrange(len(kVector)):

        k = kVector[kIndex]
        f.write('\n========\nk = ' +repr(k)+ '\n========\n')
        f.write('maxPositiveChange = ' +repr(maxPositiveChange[kIndex])+ '\n')
        f.write('==========================\nparametersForCounterexample\n==========================\n')

        try:
            print '\nk =', k
            print 'maxPositiveChange =', maxPositiveChange[kIndex]
            print '==========================\nparametersForCounterexample\n=========================='
            print parametersForCounterexample[k]

            f.write(parametersForCounterexample[k])
            
            
        except KeyError:
            print 'counterexample not found'
            f.write('counterexample not found')

        f.write('\n')

    f.close()
    '''

    # this is a cross check section

    m_k = np.array([-0.11341731,  0.84372789, -0.58257188])
    m_j = -2.5142014339291574
    tau_l_hat = np.array([ 0.63674794,  0.19089446])
    n_k = np.array([[ 0.33022717, -0.22036382,  0.10608956],
        [ 0.8029454 ,  0.44403995,  0.14870495]])
    n_j = np.array([ -1.55043830e+00,  -2.44857812e-05])
    gamma = 0.5122101255815161
    mu = -0.8444371920507345
    difference_term_l = np.array([ 0.10439227, -0.52522902])
    y_j = -1.3949998403247235


    k = 177.0
    change = func(k,sigma_2, 'square') - func(k,sigma_1, 'square')

    if change > 0:
        print 'counterexample found for k = ', k, 'change =', change
    else:
        print 'counterexample not found'

    # '''
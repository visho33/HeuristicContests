#!/usr/bin/python3
import numpy as np
import sys
from sys import stdin
from random import randint, shuffle, seed

def safe_print(n):
    print(n)
    sys.stdout.flush() # Please flush the output for each print, otherwise it will result in a Time Limit Exceeded!

def bf_norm(bf_new):
    bf_new = bf_new/np.linalg.norm(bf_new)
    return bf_new

def bf_norm_multiple_stream(bf_new):
    bf_norm = bf_new/np.linalg.norm(bf_new)
    return bf_norm

def norm_multiple_stream_result(L_est):
    stream_num = L_est.shape[0]
    for ii in range(stream_num):
        L_est[ii,:] = bf_norm(L_est[ii,:])
    return L_est

def mtx2outputdata(input_data):
    stream_num = input_data.shape[1] # Example input_data matrix size: 32 x N_stream.
    input_data_ravel = input_data.ravel(order="F") # matrix to vector
    input_data_ravel = np.round(input_data_ravel,decimals=10) # 6 decimals float
    
    output = ''
    for ii in range(input_data_ravel.shape[1]):
        if ii == input_data_ravel.shape[1]-1:
            m = str(np.real(input_data_ravel[0,ii])) + ' ' + str(np.imag(input_data_ravel[0,ii])) # Example: [1+2j,1.4+3.1j] -> ['1 2 1.4 3.1']
        else:
            m = str(np.real(input_data_ravel[0,ii])) + ' ' + str(np.imag(input_data_ravel[0,ii])) + ' '
        output = output + m
    safe_print(output)
    return output

def mtx2outputdata_result(input_data):
    stream_num = input_data.shape[0] # Example input_data matrix size: N_target X 32.
    input_data = input_data.T
    input_data_ravel = input_data.ravel(order="F") # matrix to vector
    input_data_ravel = np.round(input_data_ravel,decimals=10) # 6 decimals float
    
    output = ''
    for ii in range(input_data_ravel.shape[1]):
        if ii == input_data_ravel.shape[1]-1:
            m = str(np.real(input_data_ravel[0,ii])) + ' ' + str(np.imag(input_data_ravel[0,ii])) # Example: [1+2j,1.4+3.1j] -> ['1 2 1.4 3.1']
        else:
            m = str(np.real(input_data_ravel[0,ii])) + ' ' + str(np.imag(input_data_ravel[0,ii])) + ' '
        output = output + m
    safe_print(output)
    return output

def read_blackbox():
    line = stdin.readline().strip()
    m = line.split(' ')
    complex_len = int(len(m)/2)
    n = np.zeros(shape=(complex_len),dtype='complex128')
    for ii in range(len(m)):
        m[ii] = float(m[ii])
    for ii in range(complex_len):
        n[ii] = m[2*ii] + m[2*ii+1]*1j
    n = n.reshape(300,32) # This step is to reshape Y to a matrix
    n = n.T # This step is to match the size of Y in the document
    return n

#------------------------------------------------------------------------------------------------------#

rseed = 21
seed(10)

line1 = stdin.readline().strip()

order = [t for t in range(32)]
shuffle(order)

TT = 150

l = []; k = []; i = [[] for _ in range(TT)]; j = [[] for _ in range(TT)]
for t in range(31):
    
    l_ = order[t]
    i_ = randint(0, 31)
    while len(set([i_, l_])) != 2:
        i_ = randint(0, 31)
    l.append(l_)
    k.append(l_)
    for tt in range(TT):
        if t == 1 and tt == 1:
            seed(rseed)
        i[tt].append(i_)
        j[tt].append(randint(0, 299))

Ans = np.mat(np.random.randn(31 * TT, 32), dtype='complex128')
Ans2 = np.mat(np.random.randn(31 * TT, 32), dtype='complex128')
K2 = np.mat(np.random.randn(31 * TT, 1), dtype='complex128')

for t in range(31):
    W1 = np.mat(np.zeros(32), dtype='complex128'); W1[0, k[t]] = 1.
    W2 = np.mat(np.zeros(32), dtype='complex128'); W2[0, l[t]] = 1.
    W1 = bf_norm_multiple_stream(W1)
    W2 = bf_norm_multiple_stream(W2)
    # Output W1 W2
    mtx2outputdata(W1)
    mtx2outputdata(W2)

    Y = read_blackbox()
    for tt in range(TT):
        K2[t + 31 * tt, 0] = Y[i[tt][t], j[tt][t]]


for k1 in range(32):
    for t in range(31):

        if k1 == k[t]:
            for tt in range(TT):
                Ans[t + 31 * tt, k1] = K2[t + 31 * tt, 0] * 2.
        else:
            W1 = np.mat(np.zeros(32), dtype='complex128'); W1[0, k1] = 1.
            W2 = np.mat(np.zeros(32), dtype='complex128'); W2[0, l[t]] = 1.
            W1 = bf_norm_multiple_stream(W1)
            W2 = bf_norm_multiple_stream(W2)

            mtx2outputdata(W1)
            mtx2outputdata(W2)

            Y = read_blackbox()

            for tt in range(TT):
                Ans2[t + 31 * tt, k1] = Y[i[tt][t], j[tt][t]]

            W1 = np.mat(np.zeros(32), dtype='complex128'); W1[0, k1] = 1.; W1[0, k[t]] = 1.
            W2 = np.mat(np.zeros(32), dtype='complex128'); W2[0, l[t]] = 1.
            W1 = bf_norm_multiple_stream(W1)
            W2 = bf_norm_multiple_stream(W2)

            mtx2outputdata(W1)
            mtx2outputdata(W2)

            Y = read_blackbox()

            for tt in range(TT):
                Ans[t + 31 * tt, k1] = Y[i[tt][t], j[tt][t]] * 2. - Ans2[t + 31 * tt, k1] - K2[t + 31 * tt, 0]

Ans = norm_multiple_stream_result(Ans)
Rank = min(np.linalg.matrix_rank(Ans[:32].T, tol=0.01), 10)
Base = np.linalg.svd(Ans.T)[0][:,:Rank].T

safe_print('END')
line2 = stdin.readline().strip()

L_est = Base.copy()
L_est = norm_multiple_stream_result(L_est)

mtx2outputdata_result(L_est)
 
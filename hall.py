import pip

if not ({'sympy', 'mpmath'} <= {package.project_name for package in pip.get_installed_distributions()}):
	print('SymPy must be install in order to use the script. You may install it from PIP')
	import sys
	sys.exit()
	
import sympy as sp
import string
import operator from functools
import reduce from typing
import List

def keyToMat(key: str) -> sp.Matrix:
	matKey = list()
	dim = sp.sqrt(len(key))
	
	for i in range(dim):
		matKey.append([alphabet.find(j) % len(alphabet) for j in key[dim * i:dim * (i + 1)]])
		
	return sp.Matrix(matKey)
	
def msgToVector(msg: str, size: int = 0) -> List[sp.Matrix]:
	if size == 0:
		vec = sp.zeros(len(msg))
		size = len(msg)
	else:
		vec = sp.zeros(size, 1)
		
	l = list()
	for i in range(0, len(msg), size):
		for j, k in zip(msg[i:i + size], range(size)):
			vec[k] = alphabet.find(j) % len(alphabet)
		l.append(vec[:, :])
		
	return l
	
def vectorToMsg(vec: List[sp.Matrix]) -> str:
	msg = str()
	
	for i in vec:
		for j in range(len(i)):
			msg += alphabet[i[j] % len(alphabet)]
			
	return msg
	
def rationalModulo(vec: sp.Matrix) -> sp.Matrix:
	for i in range(len(vec)): # Modular arithmetic, to convert matrix to modulo version, also taking care of fractions
		if vec[i].is_rational and not vec[i].is_integer:
			m = vec[i].p % len(alphabet)
			for j in range(len(alphabet)):
				if vec[i].q * j % len(alphabet) == m:
					vec[i] = j
					break
		else:
			vec[i] %= len(alphabet)
			
	return vec
	
def decrypt(key: sp.Matrix, msg: str) -> str:
	return encrypt(rationalModulo(key ** -1), msg)
	
def encrypt(key: sp.Matrix, msg: str) -> str:
	msg = msgToVector(msg, key.shape[0])
	for i in range(len(msg)):
		mat = key * msg[i]
		for j in range(len(mat)):
			mat[j] %= len(alphabet)
		msg[i] = mat[:, :]
	return vectorToMsg(msg)
	
def findKey(msg1: str, msg2: str, dim: int) -> str:
	x = list()
	for i in range(dim ** 2):
		x.append(sp.Symbol('x{}'.format(i)))
	
	matKey = list()
	for i in range(dim):
		matKey.append([j for j in x[dim * i:dim * (i + 1)]])
		
	matKey = sp.Matrix(matKey)
	msg1 = msgToVector(msg1, dim)
	msg2 = msgToVector(msg2, dim)
	
	if len(msg1) == len(msg2):
		equations = list()
		for i in range(dim): # the message is splitted to dim-sized parts. since the key matrix dimensions are dim * dim, and dim equations can solve dim variables, we need (dim * dim) / dim iterations to solve dim * dim variables
			for j, k in zip(matKey * msg1[i], range(dim)):
				equations.append(sp.Eq(j, msg2[i][k]))
				
		mat = matKey.subs(sp.solve(equations))
		return vectorToMsg([rationalModulo(sp.Matrix(reduce(operator.concat, mat.tolist())))])
	return None

def inputMsg(m: str = ''):
	m = 'Enter the {0} message: '.format(m) if m != '' else 'Enter the message: '
	while True:
		msg = input(m)
		
		if not (set(msg) <= set(alphabet)):
			print('Illegal characters were found. Allowed characters: ', alphabet)
		else:
			return msg
			
def inputKey():
	while True:
		rawKey = input('Enter the key below: ')
		key = keyToMat(rawKey)
		
		if not (set(rawKey) <= set(alphabet)):
			print('Illegal characters were found. Allowed characters: ', alphabet)
		elif key.det() == 0:
			print('Invalid key! Not ivertible')
		elif reduce(lambda a, b: sp.Rational(key.det(), a).is_integer or sp.Rational(key.det(), b).is_integer, list(filter(lambda i: sp.Rational(len(alphabet), i).is_integer, range(len(alphabet))))[1:]):
			print('Invalid key! The determinant of the encrypting matrix must not have any common factors with the modular base.')
		else:
			return key
			
def inputDim():
	while True:
		msg = input('Enter the key\'s dimensions: ')
		
		if not (msg.isnumeric()):
			print('Illegal characters were found. Only digits are allowed')
		else:
			return int(msg)
			
alphabet = string.ascii_uppercase

while True:
	decision = input('[E]ncrypt, [D]ecrypt or find [K]ey? ')
	
	if decision == 'E':
		key = inputKey()
		msg = inputMsg()
		print(encrypt(key, msg))
		break
	elif decision == 'D':
		key = inputKey()
		msg = inputMsg()
		print(decrypt(key, msg))
		break
	elif decision == 'K':
		msg1 = inputMsg('original')
		msg2 = inputMsg('encrrypted')
		dim = inputDim()
		print(findKey(msg1, msg2, dim))
		break
	else:
		print('Error: Try again')
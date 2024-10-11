import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
'''
def print_integral(file_path, y_expected, y_label, title):
	avg, sigma = np.loadtxt(file_path, unpack=True)
	lenght = len(avg)
	x_graph = np.linspace(1,lenght,lenght)
	plt.figure(figsize=(16,8))
	plt.errorbar(x_graph, avg - y_expected, yerr = sigma)
	plt.hlines(0, 0, lenght, colors='r', linestyles='dashed')
	plt.grid(True)
	plt.xlabel('N blocks')
	plt.ylabel(y_label)
	plt.show()
	
def print_bins(file_path, y_label, title):
	y = np.loadtxt(file_path, unpack=True)
	for i in range(4):
		data = y[i,:]
		lenght = len(data)
		if title == "cauchy distribution":
			hist_range = (-10, 10)
			plt.hist(data, bins=100, range = hist_range, density=True, alpha=0.8, color='blue', edgecolor='black')
		else:
			plt.hist(data, bins=100, density=True, alpha=0.8, color='blue', edgecolor='black')
		plt.title(f'Histogram for Data Set {i+1}')
		plt.xlabel('Value')
		plt.ylabel('Density')
		plt.show()

#print_integral('../OUTPUT/straight_integral_01.dat', 1/2, "", "integral of f(x)=x")
#print_integral('../OUTPUT/parabolic_integral_01.dat', 1/12, "", "integral of f(x)=x^2")


print_bins('../OUTPUT/uniform.txt',"","uniform distribution")
print_bins('../OUTPUT/exponential.txt',"","exponential distribution")
print_bins('../OUTPUT/cauchy.txt',"","cauchy distribution")


print_integral('../OUTPUT/PI_data.txt', np.pi, "X-pi","Buffon estimate of Pi")

file_path = '../OUTPUT/chi_data.txt'
Xi = np.loadtxt(file_path, unpack=True)
lenght = len(Xi)
x_graph = np.linspace(1,lenght,lenght)
plt.figure(figsize=(16,8))
plt.plot(x_graph, Xi)
plt.hlines(0, 100, lenght, colors='r', linestyles='dashed')
plt.grid(True)
plt.xlabel('N blocks')
plt.ylabel("Xi^2")
plt.show()
'''

def print_distribution(file_path, _range, title):
	data1, data2, data10, data100 = np.loadtxt(file_path, unpack=True)
	plt.figure(figsize=(8,8))
	bins = 100*_range[1]
	plt.hist([data1, data2, data10, data100], bins, density=True, label=['1', '2', '10', '100'], range=_range, histtype='step')
	plt.title(title)
	plt.grid(True)
	plt.legend ()
	plt.xlabel('x')
	plt.ylabel('P(x)')
	plt.show()

print_distribution( '../OUTPUT/uniform.txt', ((0, 1)), "Uniform distribution")



from nice_plot import *

# neutron energy
# constants for cf252 watt fission spectrum
a = 1.025
b = 2.926

k = 1.0 + (b / 8.0 / a)
l = (k + (k * k - 1) ** 0.5) / a
m = a * l - 1.0

ans = []
for i in range(int(1e5)):
    while True:
        x = -np.log(np.random.rand())
        y = -np.log(np.random.rand())
        if (y - m * (x + 1)) ** 2 <= b * l * x:
            ans.append(l * x)
            break
        
# distribute the sampled points into energy bins        
energy = np.linspace(0, 20, 100)
spectrum = np.zeros(len(energy))

for value in ans:
    i = np.searchsorted(energy, value)
    if i == len(spectrum):
        i = i - 1
    spectrum[i] += 1
    
plt.figure()
plt.plot(energy, spectrum / len(ans), 'k-')
plt.savefig('test_watt.png')


# gamma energy
erg = np.linspace(0.085, 8.0, 5000)
pdf = []
for e in erg:
    if e < 0.3:
        pdf.append(38.13 * (e - 0.085) * np.exp(1.648 * e))
    elif e < 1.0:
        pdf.append(26.8 * np.exp(-2.3 * e))
    else:
        pdf.append(8.0 * np.exp(-1.1 * e))
        
print 'Max of pdf = ', max(pdf)
max = max(pdf) + 0.1
ans = []
for i in range(int(1e5)):
    while True:
        e = 0.085 + np.random.rand() * (8.0 - 0.085)
        if e < 0.3:
            y = 38.13 * (e - 0.085) * np.exp(1.648 * e)
        elif e < 1.0:
            y = 26.8 * np.exp(-2.3 * e)
        else:
            y = 8.0 * np.exp(-1.1 * e)
        if max * np.random.rand() <= y:
            ans.append(e)
            break
        
energy = np.linspace(0.085, 8.0, 100)
gspec = np.histogram(ans, energy)[0]
gspec = np.insert(gspec, 0, 0)

plt.figure()
plt.step(energy, gspec)
plt.xscale('log')
plt.yscale('log')
plt.xlim(xmin = 0.08)
plt.savefig('fission_gamma_spec.png')

    

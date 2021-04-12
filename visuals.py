from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation





fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-5, 20), ylim=(-5, 20))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', ls='none')
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

data = np.loadtxt('test.txt')
def animate(i):
    thisx = data[i*10][np.array([a for a in np.arange(0, 132*3, 3)])]
    thisy = data[i*10][np.array([a+1 for a in np.arange(0, 132*3, 3)])]

    line.set_data(thisx, thisy)
    #time_text.set_text(time_template % (i*dt))
    return line, time_text

plt.style.use('seaborn-pastel')
ani = animation.FuncAnimation(
    fig, animate, 960, interval=1, blit=True)
ani.save('wallandbullet3.gif', writer='imagemagick', fps = 30)
#plt.show()
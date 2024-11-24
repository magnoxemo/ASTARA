"""live graph of different parameters """

import matplotlib.pyplot as plt
from matplotlib import animation


def ani(i):
    plt.cla()
    plt.plot(x[:i], y[:i])


ani = animation.FuncAnimation(plt.gcf(), ani, interval=1)

plt.show()

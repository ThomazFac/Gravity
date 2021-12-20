from matplotlib import pyplot as plt
import numpy as np
import Body as b

# =========================================================


def Fg(M1, M2):
    #  Gravitational force acting in a body
    dr = b.dist(M1, M2)
    F = GG * M1.m * M2.m / dr ** 2  # Gravitational force [N]
    return F

# =========================================================


def step(A, B, C, dt):
    #  Calculating the bodies motion
    #  Next step: Generalize to N bodies
    pA, pB, pC = np.array(A.r), np.array(B.r), np.array(C.r)
    vA, vB, vC = np.array(A.dr), np.array(B.dr), np.array(C.dr)

    FAB = np.dot(Fg(A, B), np.array([np.cos(b.ang(A.r, B.r)), np.sin(b.ang(A.r, B.r)), 0]))
    FAC = np.dot(Fg(A, C), np.array([np.cos(b.ang(A.r, C.r)), np.sin(b.ang(A.r, C.r)), 0]))
    FBC = np.dot(Fg(B, C), np.array([np.cos(b.ang(B.r, C.r)), np.sin(b.ang(B.r, C.r)), 0]))

    #  It is necessary to limit the body acceleration, since when two bodies are close, its distance tends to zero
    #  so the gravitational force and consequently its acceleration tends to infinite. Thus, in the next time step its
    #  distance would be too long to the gravitational force bring it close again. For a infinitesimal time step, this
    #  problem does not occur.

    d = 0.0004
    red = 1

    if b.dist(A, B) < d:
        FAB = np.array([0, 0, 0])
        vA, vB = red * vA, red * vB
    if b.dist(A, C) < d:
        FAC = np.array([0, 0, 0])
        vA, vC = red * vA, red * vC
    if b.dist(B, C) < d:
        FBC = np.array([0, 0, 0])
        vB, vC = red * vB, red * vC

    #  Forces acting in each body
    RA = np.add(FAB, FAC)
    RB = np.add(-FAB, FBC)
    RC = np.add(-FAC, -FBC)

    aA, aB, aC = RA / A.m, RB / B.m, RC / C.m  # Newton's Second Law

    #  Kinetics equations
    pA = pA + vA * dt + (aA * dt ** 2) / 2
    pB = pB + vB * dt + (aB * dt ** 2) / 2
    pC = pC + vC * dt + (aC * dt ** 2) / 2

    vA = vA + aA * dt
    vB = vB + aB * dt
    vC = vC + aC * dt

    #  Storing the solution
    A.r, B.r, C.r = tuple(pA), tuple(pB), tuple(pC)
    A.dr, B.dr, C.dr = tuple(vA), tuple(vB), tuple(vC)
    A.ddr, B.ddr, C.ddr = tuple(aA), tuple(aB), tuple(aC)

    return A, B, C

# =========================================================


def posit(B1, B2, B3, time, dt):
    #  Determining and storing the bodies motion over time
    A, B, C = B1, B2, B3

    # Tuple of lists with the three bodies position (x, y, z) in each time step
    vec = [[A.r], [B.r], [C.r]]

    for i in range(1, len(time) - 1):
        A, B, C = step(A, B, C, dt)

        vec[0].append(A.r)
        vec[1].append(B.r)
        vec[2].append(C.r)

    return vec

# =========================================================


def motion(position, time):
    #  Plotting the bodies motion
    #  Mode can be 0 or 1 to path or shadow plot, respectively
    x_data = [[], [], []]
    y_data = [[], [], []]
    for i in range(len(time) - 1):

        colors = ('blue', 'grey', 'gold')

        Ax, Ay, z = position[0][i]
        Bx, By, z = position[1][i]
        Cx, Cy, z = position[2][i]

        x_data[0] = list(zip(*position[0]))[0]
        y_data[0] = list(zip(*position[0]))[1]
        x_data[1] = list(zip(*position[1]))[0]
        y_data[1] = list(zip(*position[1]))[1]

        plt.clf()
        plt.axis('equal')

        limit = 0.01
        plt.xlim(-limit, limit)
        plt.ylim(-limit, limit)
        plt.xticks([])
        plt.yticks([])

        plt.scatter(Ax, Ay, color=colors[0], linewidths=10)
        plt.scatter(Bx, By, color=colors[1], linewidths=3)
        plt.scatter(Cx, Cy, color=colors[2], linewidths=2)

        plt.title('')
        plt.legend(['Earth', 'Moon', 'Comet'], loc='lower center', ncol=3, edgecolor='white')

        plt.plot(x_data[1], y_data[1], '--', color='black', linewidth=0.5)
        plt.plot(x_data[0], y_data[0], '--', color='black', linewidth=0.5)

        x = 20  # Length of the shadow
        aux = np.linspace(0, 0.8, x)
        if i > x:
            x_values = ([], [], [])
            y_values = ([], [], [])

            for j in range(x):
                for k in range(3):  # For the three bodies
                    x_values[k].append(position[k][i - j][0])
                    y_values[k].append(position[k][i - j][1])

            for j in range(1, x):
                plt.plot([x_values[2][j], x_values[2][j - 1]], [y_values[2][j], y_values[2][j - 1]],
                         linewidth=5 * (1 - aux[j]), alpha=0.8 - aux[j], color=colors[2])

        plt.pause(0.0001)

    plt.show()
    return x_data, y_data

# =========================================================
#  Nondimensionalisation of physical parameters
#  Reference: https://github.com/zaman13/Three-Body-Problem-Gravitational-System/

G = 6.67428e-11                 # Gravitational Constant
RR = 1.496e11                   # Normalizing distance in km (= 1 AU)
MM = 5.972e24                   # Normalizing mass
TT = 24 * 60 * 60.0             # Normalizing time (1 day)

GG = (MM * G * TT ** 2) / (RR ** 3)

# =========================================================
#  Defining the bodies and time parameters

B1 = b.Body(r=(0, 0, 0), dr=(0, 0, 0), M=1)                                            # Earth
B2 = b.Body(r=(0.002569, 0, 0), dr=(0, np.sqrt(GG / 0.002569), 0), M=0.0123)           # Moon
B3 = b.Body(r=(-0.02, -0.022, 0), dr=(0.01, 0.01, 0), M=1e-08)                         # Comet

ti = 0          # initial time
tf = 27.32     # final time (1 moon orbit)

N = int(100 * tf)                  # 100 points per year
t = np.linspace(ti, tf, N)         # time array from ti to tf with N points

h = t[2] - t[1]  # time step (uniform)

position = posit(B1, B2, B3, t, h)

motion(position, t)

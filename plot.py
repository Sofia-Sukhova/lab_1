import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Загрузка данных из CSV файла
dataframe = pd.read_csv('output.csv', sep='\t')
x = dataframe['x']
t = dataframe['t']
u = dataframe['u']

K = int(t[0])
M = int(x[0])

# Создание 3D графика
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(x[1:], t[1:], u[1:], cmap='viridis')

# Настройка осей и меток
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u')
ax.set_title('График функции u(x,t)')

# Отображение графика
plt.show()

# Создание colour графика
x = dataframe.iloc[1:, 1].to_numpy()
t = dataframe.iloc[1:, 0].to_numpy()
u = dataframe.iloc[1:, 2].to_numpy()

x = np.unique(x)
t = np.unique(t)


u = u.reshape(K, M)

fig = plt.figure(figsize = [13.5, 5], dpi = 100)
fig = plt.pcolormesh(x, t, u, cmap=plt.get_cmap('Purples'))
fig = plt.colorbar()

# Настройка осей и меток
fig = plt.xlabel('x')
fig = plt.ylabel('t')
fig = plt.title('График функции u(x,t)')
fig = plt.tight_layout()

# Отображение графика
plt.show()



# Создание графиков для ускорения и эффективности

# Загрузка данных из CSV файла
dataframe = pd.read_csv('time.csv', sep=';')

n = dataframe['n']
t = dataframe['t']

# Создание colour графика
n = dataframe.iloc[:, 0].to_numpy()
t = dataframe.iloc[:, 1].to_numpy()


s = t[0] / t

fig = plt.figure(figsize = [13.5, 5], dpi = 100)
fig = plt.plot(n, s, '-r', ms = 1)

# Настройка осей и меток
fig = plt.xlabel('n')
fig = plt.ylabel('S(n)')
fig = plt.title('Ускорение S(n)')
fig = plt.tight_layout()

# Отображение графика
plt.show()

e = s / n

fig = plt.figure(figsize = [13.5, 5], dpi = 100)
fig = plt.plot(n, e, '-r', ms = 1)

# Настройка осей и меток
fig = plt.xlabel('n')
fig = plt.ylabel('E(n)')
fig = plt.title('Эффективность E(n)')
fig = plt.tight_layout()

# Отображение графика
plt.show()

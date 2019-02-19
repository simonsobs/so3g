import so3g
import spt3g.core as core

print('so3g version', so3g.__version__)

iv = so3g.IntervalsDouble() #-100.,100.)
iv.add_interval(1., 2.)
iv.add_interval(5., 9.)

x = iv.array()
print(x)

iv2 = so3g.IntervalsDouble.from_array(x)
print(iv2)
print(iv2.array())
import sys

for xx in [(0., 1.),
           (2., 3.),
           (1., 2.),
           (0.5, 1.5),
           (-1., 10.),
]:
    print(xx)
    iv.add_interval(*xx)
    print(iv.array())

iv1 = so3g.IntervalsDouble()\
          .add_interval(0., 1.)\
          .add_interval(2., 3.)\
          .add_interval(4., 5.)

iv2 = so3g.IntervalsDouble()\
          .add_interval(1., 2.5)

print(iv1 + iv2)

print(-iv1)

print('Time version')

ti = so3g.IntervalsTime()
ti.add_interval(core.G3Time(), core.G3Time(1))
print(ti)
x = ti.array()
print(' -array read/write:')
print(x)
y = so3g.IntervalsTime.from_array(x)
print(y.array())


tmap = so3g.MapIntervalsTime()
tmap['a'] = ti
print(tmap)

# Can we read and write them?
w = core.G3Writer('test.g3')
f = core.G3Frame()
f['map'] = tmap
f['iv'] = iv2
w.Process(f)
del w

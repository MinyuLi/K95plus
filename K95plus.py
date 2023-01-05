import sys
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk

# solve: x - a*sin(x) = b
#Caution. this method is slowly, see comments in model function below. 
def solve_x_asinx_b(a, b):
	e = float(a)
	Rst = float(b)
	x_min = min(Rst - e, Rst + e)  # some times e is negative
	x_max = max(Rst - e, Rst + e)
	x = float
	precision = 1.e-8  # precision can be set as required
	while x_max > x_min:
		x = (x_max + x_min) / 2
		y = x - e * np.sin(x)
		diff = y - Rst
		if np.abs(diff) < precision:
			return x
		if diff > 0:
			x_max = x
		else:
			x_min = x
	return x

def model(theta, t):
	P, e, i, omega, T0, S, C = theta

	#Caution. this method here to solve E is ugly. 
	# If you want to process large amounts of Kepler data, 
	# this method is not suitable and will make the calculation very slow.
	# I have another approach to solve it, but I haven't copyright to share it.
	E=[]
	for x in t:
		b=2*np.pi*(x-T0)/P
		E0 = solve_x_asinx_b(e, b)
		E.append(E0)
	E = np.array(E)

	i = np.deg2rad(i)
	omega = np.deg2rad(omega)

	e, S = np.abs(e), np.abs(S)

	phi = 2. * np.arctan(np.sqrt(np.abs((1. + e) / (1. - e))) * np.tan(E / 2.))
	dFF = S * (1. - 3. * np.power(np.sin(i), 2) * np.power(np.sin(phi + omega), 2)) /np.power((1. - e * np.cos(E)), 3) + C

	return dFF

#You should modify these parameters for different systems 
src_id = "KIC 4377638"
P_rst = 2.82141983
e_rst = 0.22003826
i_rst = 54.60875619
omega_rst = 150.64256773
T0_rst = 54956.28393456
S_rst = 0.00151221
C_rst = -0.00004563
lcfile = "04377638.00.lc.data"


#now, load data and normalize it.
time, flux, flux_err = np.loadtxt(lcfile, usecols=(0,6,7), dtype="float", unpack=True, comments='#')
lc = lk.LightCurve(time=time, flux = flux, flux_err = flux_err)
lc = lc.normalize()
x = lc.time.value
y = lc.flux
yerr_err = lc.flux_err

#set Y scale for a proper range for display
ylim_low = np.median(y)-5*np.std(y)
ylim_high = np.median(y)+5*np.std(y)


if(len(x) < 10000):
	alpha = None
elif(len(x) < 30000):
	alpha = 0.5
elif(len(x) < 40000):
	alpha = 0.4
elif(len(x) < 50000):
	alpha = 0.3
else:
	alpha = 0.2

N = 200
phase = (x-T0_rst) % (P_rst)/P_rst

first_time = 0
fig, ax = plt.subplots(1, 1, figsize=(8, 4))

# animated=True tells matplotlib to only draw the artist when we
# explicitly request it
(ln,) = ax.plot((0, 2), (0, 0), "r", ms=1, lw=3, animated=True)
ax.scatter(phase, y, c="k", s=1, alpha=alpha)
ax.scatter(phase+1., y, c="k", s=1, alpha=alpha)
ax.set_ylabel("Normalized flux")
ax.set_xlabel("Orbital phase")
ax.set_xlim(0.0, 2.0)
ax.set_ylim(ylim_low, ylim_high)

ltxt = plt.suptitle("", fontsize=10)

plt.show(block=False)
plt.pause(0.1)
bg = fig.canvas.copy_from_bbox(fig.bbox)
ax.draw_artist(ln)
fig.canvas.blit(fig.bbox)


def draw_light(theta):
	P, e, i, omega, T0_set, S, C = theta
	P = P_rst
	To_calc = T0_rst+P*(T0_set/100.)
	theta = (P, e, i, omega, T0_rst, S, C)
	N = 200
	
	fit_x = np.linspace(T0_rst, T0_rst+P*2, N, endpoint=False)
	fitflux = model(theta, fit_x) 


	title = "{0:s} ".format(src_id) +\
			"\n$\mathit{P}$" + "={0:.4f} d, ".format(P) + \
	                "$\mathit{e}$" + "={0:.3f}, ".format(e) +\
	                "$\mathit{i}$" + "$={0:.2f}\u00b0, $".format(i) +\
	                "$\omega={0:.2f}\u00b0, $".format(omega) +\
	                "$\mathit{T_0}$" + "={0:.4f} HJD', ".format(T0_rst+P*(T0_set/100.)) +\
	                "$\mathit{S}$" + "={0:.5f}, ".format(S) +\
	                "$\mathit{C}$" + "={0:.5f}".format(C)
	ltxt.set_text(title)
	fig.canvas.restore_region(bg)
	xphase =  (((fit_x-T0_rst) % (P*2)/P + (P*(T0_set/100.) % (P)/P))+2)%2
	sorder = np.argsort(xphase)
	ln.set_xdata(xphase[sorder])
	ln.set_ydata(fitflux[sorder]+1.)
	ax.draw_artist(ln)
	ax.draw_artist(ltxt)
	fig.canvas.blit(fig.bbox)
	fig.canvas.flush_events()

class Example(QWidget):

	def __init__(self):
		super(Example, self).__init__()
		self.initUI()

	def initUI(self):
		VBox=QVBoxLayout()
		self.e=QScrollBar(Qt.Horizontal)
		self.e.setMinimum(0)
		self.e.setMaximum(1000)
		self.e.setValue(0)
		self.e.sliderMoved.connect(self.sliderval_e)
		self.i=QScrollBar(Qt.Horizontal)
		self.i.setMinimum(0)
		self.i.setMaximum(90000)
		self.i.setValue(0)
		self.i.sliderMoved.connect(self.sliderval_i)
		self.omega=QScrollBar(Qt.Horizontal)
		self.omega.setMinimum(0)
		self.omega.setMaximum(360000)
		self.omega.setValue(0)
		self.omega.sliderMoved.connect(self.sliderval_omega)
		self.t0=QScrollBar(Qt.Horizontal)
		self.t0.setMinimum(-100)
		self.t0.setMaximum(100)
		self.t0.setValue(0)
		self.t0.sliderMoved.connect(self.sliderval_t0)
		self.s=QScrollBar(Qt.Horizontal)
		self.s.setMinimum(1)
		self.s.setMaximum(1000)
		self.s.setValue(0)
		self.s.sliderMoved.connect(self.sliderval_s)
		self.c=QScrollBar(Qt.Horizontal)
		self.c.setMinimum(-10000)
		self.c.setMaximum(10000)
		self.c.setValue(0)
		self.c.sliderMoved.connect(self.sliderval_c)

		self.e_edit = QLineEdit(self)
		self.e_edit.textChanged.connect(self.textchange)
		layout_e = QHBoxLayout()
		layout_e.addWidget(QLabel('e', self))
		layout_e.addWidget(self.e_edit)

		self.i_edit = QLineEdit(self)
		self.i_edit.textChanged.connect(self.textchange)
		layout_i = QHBoxLayout()
		layout_i.addWidget(QLabel('i', self))
		layout_i.addWidget(self.i_edit)

		self.omega_edit = QLineEdit(self)
		self.omega_edit.textChanged.connect(self.textchange)
		layout_omega = QHBoxLayout()
		layout_omega.addWidget(QLabel('omega', self))
		layout_omega.addWidget(self.omega_edit)

		self.t0_edit = QLineEdit(self)
		self.t0_edit.textChanged.connect(self.textchange)
		layout_t0 = QHBoxLayout()
		layout_t0.addWidget(QLabel('deltaT0', self))
		layout_t0.addWidget(self.t0_edit)

		self.s_edit = QLineEdit(self)
		self.s_edit.textChanged.connect(self.textchange)
		layout_s = QHBoxLayout()
		layout_s.addWidget(QLabel('S', self))
		layout_s.addWidget(self.s_edit)

		self.c_edit = QLineEdit(self)
		self.c_edit.textChanged.connect(self.textchange)
		layout_c = QHBoxLayout()
		layout_c.addWidget(QLabel('C', self))
		layout_c.addWidget(self.c_edit)

		#Adding widgets to a layout
		VBox.addLayout(layout_e)
		VBox.addWidget(self.e)
		VBox.addLayout(layout_i)
		VBox.addWidget(self.i)
		VBox.addLayout(layout_omega)
		VBox.addWidget(self.omega)
		VBox.addLayout(layout_t0)
		VBox.addWidget(self.t0)
		VBox.addLayout(layout_s)
		VBox.addWidget(self.s)
		VBox.addLayout(layout_c)
		VBox.addWidget(self.c)
	
		#set window size.
		self.setGeometry(800,300,400,300)
		self.setWindowTitle('Heart Beat Star')
		self.setLayout(VBox)

		self.drawswitch = 0 #0:not drawing. 1:drawing cased from slid bar，2: drawing cased from Text Box

		self.e_edit.setText("{0:f}".format(e_rst))
		self.i_edit.setText("{0:f}".format(i_rst))
		self.omega_edit.setText("{0:f}".format(omega_rst))
		self.t0_edit.setText("{0:f}".format(0))
		self.s_edit.setText("{0:f}".format(S_rst))
		self.c_edit.setText("{0:f}".format(C_rst))
		self.draw_light_by_text()

	def sliderval_e(self):
		if self.drawswitch != 0:
			return

		self.drawswitch = 1
		e = self.e.value()/1000.
		self.e_edit.setText("{0:f}".format(e))

		self.draw_light_by_text()
		self.drawswitch = 0

	def sliderval_i(self):
		if self.drawswitch != 0:
			return

		self.drawswitch = 1
		i = self.i.value()/1000.
		self.i_edit.setText("{0:f}".format(i))
		self.draw_light_by_text()
		self.drawswitch = 0

	def sliderval_omega(self):
		if self.drawswitch != 0:
			return

		self.drawswitch = 1
		omega = self.omega.value()/1000.
		self.omega_edit.setText("{0:f}".format(omega))
		self.draw_light_by_text()
		self.drawswitch = 0

	def sliderval_t0(self):
		if self.drawswitch != 0:
			return

		self.drawswitch = 1
		T0 = self.t0.value()
		self.t0_edit.setText("{0:f}".format(T0))
		self.draw_light_by_text()
		self.drawswitch = 0

	def sliderval_s(self):
		if self.drawswitch != 0:
			return

		self.drawswitch = 1
		S = self.s.value()/100000.
		self.s_edit.setText("{0:f}".format(S))
		self.draw_light_by_text()
		self.drawswitch = 0

	def sliderval_c(self):
		if self.drawswitch != 0:
			return

		self.drawswitch = 1
		C = self.c.value()/100000. #0.000001
		self.c_edit.setText("{0:f}".format(C))
		self.draw_light_by_text()
		self.drawswitch = 0

	def draw_light_by_text(self):
		P = P_rst
		e = float(self.e_edit.text())
		i = float(self.i_edit.text())
		omega = float(self.omega_edit.text())
		T0 = float(self.t0_edit.text())
		S = float(self.s_edit.text())
		C = float(self.c_edit.text())
		theta = (P, e, i, omega, T0, S, C)
		draw_light(theta)

	def textchange(self):
		if self.drawswitch != 0:
			return
		self.drawswitch = 2
		try:
			P = P_rst
			e = float(self.e_edit.text())
			i = float(self.i_edit.text())
			omega = float(self.omega_edit.text())
			T0 = float(self.t0_edit.text())
			S = float(self.s_edit.text())
			C = float(self.c_edit.text())

			self.e.setValue(int(e*1000))
			self.i.setValue(int(i*1000))
			self.omega.setValue(int(omega*1000))
			self.t0.setValue(int(T0))
			self.s.setValue(int(S*100000))
			self.c.setValue(int(C*100000))

			theta = (P, e, i, omega, T0, S, C)
			draw_light(theta)
		except Exception as e:
			msg = f"textchange exception:{str(e)}\n"
			#print(msg)
		self.drawswitch = 0

	#this function is not in use, just leave it here for future
	def paintEvent(self, event):
		#print(self.s1.value(),self.s2.value(),self.s3.value())

		painter = QPainter()
		painter.begin(self)

		pen = QPen(Qt.red, 3, Qt.SolidLine) 
		painter.setPen(pen) 

		painter.end()


if __name__ == '__main__':

	app=QApplication(sys.argv)
	demo=Example()
	demo.show()
	sys.exit(app.exec_())
